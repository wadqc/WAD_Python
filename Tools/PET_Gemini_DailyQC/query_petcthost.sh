#!/bin/bash
# script for querying the Philips Gemini PET-CT HOST computer for daily QC data
# to be used inside cron (crontab -e), e.g.:
#
#   # check every half an hour between 8:00-17:00 on working days for new QC data
#   0,30 08-17 * * 1-5 /etc/scripts/query_petcthost.sh 2>&1 > /dev/null
#
# requirements: dcmtk installed under Linux

SOURCE_IP=petct-host
SOURCE_PORT=104
SOURCE_AET=PETCT-HOST

TARGET_IP=wadserver
TARGET_PORT=11112
TARGET_AET=DCM4CHEE

PID="r057"
NAME="PET*"
CONTENTDATE=$(date +"%Y%m%d")
# prevent retrieving 6:00 background measurements
CONTENTTIME="063000-235900"
SERIESDESCRIPTION="QC*"


QUERY=$(findscu -S -aet $TARGET_AET -aec $SOURCE_AET -k 0008,0052="STUDY" -k PatientID="$PID" -k PatientName="$NAME" -k StudyInstanceUID $SOURCE_IP $SOURCE_PORT 2>&1)

if [ $? -ne 0 ]
then
   echo failed!
   continue
fi

NRMATCHES=$(echo "$QUERY" | grep Response | wc -l)

if [ $NRMATCHES -ne 1 ]
then
   echo
   echo No studies found!
   continue
else
   echo
   StudyInstanceUID=$(echo "$QUERY" | grep StudyInstanceUID | cut -f1 -d] | cut -f2 -d[)
   echo StudyInstanceUID = $StudyInstanceUID
   QUERY=$(findscu -S -aet $TARGET_AET -aec $SOURCE_AET -k 0008,0052="SERIES" -k StudyInstanceUID="$StudyInstanceUID" -k SeriesInstanceUID -k ContentDate="$CONTENTDATE" -k ContentTime="$CONTENTTIME" -k SeriesDescription="$SERIESDESCRIPTION" $SOURCE_IP $SOURCE_PORT 2>&1)

   if [ $? -ne 0 ]
   then
      echo failed!
      continue
   fi

   NRMATCHES=$(echo "$QUERY" | grep Response | wc -l)
   echo $NRMATCHES series found on source
   SeriesInstanceUIDs_SOURCE=$(echo "$QUERY" | grep SeriesInstanceUID | cut -f1 -d] | cut -f2 -d[)
#   for i in $SeriesInstanceUIDs_SOURCE
#   do
#      echo Serie on source: $i
#   done

fi





QUERY=$(findscu -S -aet $TARGET_AET -aec $TARGET_AET -k 0008,0052="SERIES" -k StudyInstanceUID="$StudyInstanceUID" -k SeriesInstanceUID -k ContentDate="$CONTENTDATE" -k ContentTime="$CONTENTTIME" -k SeriesDescription="$SERIESDESCRIPTION" $TARGET_IP $TARGET_PORT 2>&1)

   if [ $? -ne 0 ]
   then
      echo failed!
      continue
   fi

   SeriesInstanceUIDs_TARGET=$(echo "$QUERY" | grep SeriesInstanceUID | cut -f1 -d] | cut -f2 -d[)

   for i in $SeriesInstanceUIDs_SOURCE
   do
      match=$(echo "$SeriesInstanceUIDs_TARGET" | grep $i | wc -l)
      if [ $match -ne 1 ]
      then
         echo $i not yet present on target
         movescu -v -S -aem $TARGET_AET -aec $SOURCE_AET -aet $TARGET_AET +xa -k 0008,0052="SERIES" -k PatientID="$PID" -k StudyInstanceUID="$StudyInstanceUID" -k SeriesInstanceUID="$i" $SOURCE_IP $SOURCE_PORT
      else
         echo $i present on target
      fi
   done
