#!/usr/bin/env python
"""
http://www.protocolostomy.com/2008/10/22/generating-reports-with-charts-using-python-reportlab/
http://www.reportlab.com/software/documentation/

TODO:
    * kleuren aanpassen (lichtgeel onzichtbaar)
    * critmin/max dikke lijnen erin
"""
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import pagesizes
from reportlab.lib.units import inch
from reportlab.graphics.charts.legends import Legend
from reportlab.graphics.samples.excelcolors import *
from reportlab.lib import colors

import time

from reportlab.pdfgen import canvas
from reportlab.lib.units import mm

class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        canvas.Canvas.__init__(self, *args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        """add page info to each page (page x of y)"""
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_number(num_pages)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        # Change the position of this to wherever you want the page number to be
        self.drawRightString(185*mm, 72.0-5 * mm ,
                             "%d/%d" % (self._pageNumber, page_count))


def dateformatter(val):
    #return val
    #return 'x=%s'%val
    return time.strftime('%Y-%m-%d', time.localtime(val))

def emptyformatter(val):
    return ''

class Reporter:
    defaultPageSize = None
    PAGE_HEIGHT = None
    styles = None
    Elements = None
    HeaderStyle  = None
    ParaStyle = None
    PreStyle = None
    colors = []
    maxnumgraphsinchart = 9
    def defaultstyle(self):
        self.defaultPageSize = pagesizes.A4
        self.PAGE_HEIGHT=self.defaultPageSize[1]
        self.PAGE_WIDTH = self.defaultPageSize[0]
        self.styles = getSampleStyleSheet()
        self.Elements=[]
        self.HeaderStyle = self.styles["Heading1"]
        self.Header2Style = self.styles["Heading2"]
        self.ParaStyle = self.styles["Normal"]
        self.PreStyle = self.styles["Code"]
#        self.colors = 5*[color01,color02,color03,color04,color05,color06,color07,color08,color09,color10,
#                       color01Light,color02Light,color03Light,color04Light,color05Light,color06Light,color07Light,color08Light,color09Light,color10Light,
#                       color01Dark,color02Dark,color03Dark,color04Dark,color05Dark,color06Dark,color07Dark,color08Dark,color09Dark,color10Dark]
        self.colorRangeStatic(64)

    def __init__(self,title='Sample report',author='Q.C. Wad', url='http://localhost:8080',abstract='This is a sample report'):
        self.defaultstyle()
        self.Title = title
        self.Author = author
        self.URL = url
        self.Abstract = abstract

    @staticmethod
    def _header_footer(canvas, doc):
        # Save the state of our canvas so we can draw on it
        canvas.saveState()
        styles = getSampleStyleSheet()

        # Header
        ##header = Paragraph('This is a multi-line header.  It goes on every page.   ' * 5, styles['Normal'])
        ##w, h = header.wrap(doc.width, doc.topMargin)
        ##header.drawOn(canvas, doc.leftMargin, doc.height + doc.topMargin - h)

        # Footer
        canvas.line(doc.leftMargin,doc.bottomMargin,doc.leftMargin+doc.width,doc.bottomMargin)
        footer = Paragraph('QC Status Rapportage', styles['Italic'])
        w, h = footer.wrap(doc.width, doc.bottomMargin-5*mm)
        footer.drawOn(canvas, doc.leftMargin, doc.bottomMargin-5*mm)

        # Release the canvas
        canvas.restoreState()

    def header(self,txt, style=None, klass=Paragraph, sep=0.3):
        if style is None:
            style = self.HeaderStyle

        s = Spacer(0.2*inch, sep*inch)
        para = klass(txt, style)
        sect = [s, para]
        result = KeepTogether(sect)
        return result

    def paragraph(self,txt):
        return self.header(txt, style=self.ParaStyle, sep=0.1)

    def header2(self,txt):
        return self.header(txt, style=self.Header2Style, sep=0.1)

    def preblock(self,txt):
        s = Spacer(0.1*inch, 0.1*inch)
        p = Preformatted(txt, self.PreStyle)
        precomps = [s,p]
        result = KeepTogether(precomps)
        return result

    def render(self,outputname='gfe.pdf'):
        doc = SimpleDocTemplate(outputname)
        doc.build(self.Elements,onFirstPage=self._header_footer, onLaterPages=self._header_footer,
                  canvasmaker=NumberedCanvas)

    def addtitlepage(self):
        self.Elements.insert(0,self.header(self.Title))
        self.Elements.insert(1,self.paragraph(self.Author))
        self.Elements.insert(2,self.paragraph(self.URL))
        self.Elements.insert(3,self.header2("Samenvatting"))
        self.Elements.insert(4,self.paragraph(self.Abstract))

    def addtitlepagewithtable(self,tabledata):
        self.Elements.insert(0,self.header(self.Title))
        self.Elements.insert(1,self.paragraph(self.Author))
        self.Elements.insert(2,self.paragraph(self.URL))
        self.Elements.insert(3,self.header2("Samenvatting"))
        self.Elements.insert(4,self.table(tabledata,lines=False,align='LEFT'))

    def addsection(self,sectionlist):
        self.Elements.append(KeepTogether(sectionlist))


    def table(self,data,lines=True,align='CENTER'):
        s = Spacer(0.1*inch, 0.1*inch)
        dtable = Table(data)
        dtable.hAlign= align
        if lines:
            dtable.setStyle(TableStyle([('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
                                        ('BOX', (0, 0), (-1, -1), 0.25, colors.black)]))
        result = (KeepTogether([s,dtable]))
        return result

    def dategraph(self,lxydata,tmin=None,tmax=None): # list of (label, [x,yarr]) with x in seconds
        maxX = lxydata[0][1][0][0]
        maxY = 1.
        minX = maxX
        minY = -1
        for label,xyarr in lxydata:
            xarr = [x for (x,y) in xyarr]
            yarr = [y for (x,y) in xyarr]
            maxY = max(maxY, max(yarr))
            minY = min(minY, min(yarr))
            maxX = max(maxX, max(xarr))
            minX = min(minX, min(xarr))
        if not tmin is None:
            minX = tmin
        if not tmax is None:
            maxX = tmax
        lxydata.append(('CritMax',[[minX,1.],[maxX,1.]]))
        lxydata.append(('CritMin',[[minX,-1.],[maxX,-1.]]))

        drawing = Drawing(self.PAGE_WIDTH-150, 200)
        lc = LinePlot()
        lc.x = 30
        lc.y = 50
        lc.height = 125
        lc.width = self.PAGE_WIDTH-150-30
        lc.data = [xy for (l,xy) in lxydata]
        datalabels = [ l for (l,xy) in lxydata]

        lc.xValueAxis.valueMin = minX
        lc.xValueAxis.valueMax = maxX
        lc.xValueAxis.valueStep = int((maxX-minX)/5)
        lc.xValueAxis.labelTextFormat = dateformatter # Use the formatter
        lc.xValueAxis.visibleGrid           = True
        lc.xValueAxis.gridStrokeWidth = .5
        lc.xValueAxis.gridStrokeColor = colors.darkgrey
        lc.xValueAxis.subTickNum = 1
        lc.xValueAxis.subGridStrokeWidth = .25
        lc.xValueAxis.subGridStrokeColor = colors.lightgrey
        lc.xValueAxis.visibleSubGrid = True

        lc.yValueAxis.valueMin = minY
        lc.yValueAxis.valueMax = maxY
        lc.yValueAxis.valueStep = .25
        if (maxY-minY)>2:
            lc.yValueAxis.valueStep = .5
        if (maxY-minY)>4:
            lc.yValueAxis.valueStep = 1.
        lc.yValueAxis.visibleGrid           = True
        lc.yValueAxis.gridStrokeWidth = .5
        lc.yValueAxis.gridStrokeColor = colors.darkgrey
        lc.yValueAxis.subTickNum = 1
        lc.yValueAxis.subGridStrokeWidth = .25
        lc.yValueAxis.subGridStrokeColor = colors.lightgrey
        lc.yValueAxis.visibleSubGrid = True

        for x in range(len(lc.data)):
            lc.lines[x].strokeColor = self.colors[x]
            lc.lines[x].name = datalabels[x]
            if len(lc.data[x]) ==1:
                lc.lines[x].symbol = makeMarker('FilledCircle') # added to make filled circles.
            lc.lines[x].strokeWidth = 2
            if lc.lines[x].name == 'CritMin' or lc.lines[x].name == 'CritMax': # distinguish min/max
                lc.lines[x].strokeColor = self.colors[0]
                lc.lines[x].strokeDashArray = (5, 1)

        self.addLabels(drawing,ylabel="geschaalde waarde")
        drawing.add(lc)
        drawing.add(self.addAutoLegend(lc,(lc.x,lc.y),len(lc.data)))
        return drawing

    def addAutoLegend(self,chart,offset,num,side='bottom'):
        from reportlab.lib.validators import Auto
        width = 300
        height = 150
        legend = Legend()
        if side == 'bottom':
            legend.x = offset[0]-15
            legend.dx = 20
            legend.dy = 5
            legend.y = offset[1]-25
            legend.deltax = None
            legend.autoXPadding = 35
            legend.deltay = 5
            legend.boxAnchor = 'nw'
            legend.dxTextSpace = 5
            legend.columnMaximum = (num/3)
            if num%3 >0:
                legend.columnMaximum += 1

            legend.variColumn = True
        elif side == 'right':
            legend.x = offset[0]+offset[2]+40
            legend.dx = 20
            legend.dy = 5
            legend.y = offset[1]+offset[3]-15
            legend.deltax = None
            legend.deltay = 5
            legend.boxAnchor = 'nw'
            legend.dxTextSpace = 5
            legend.columnMaximum = 9

        legend.colorNamePairs=Auto(chart=chart)
        return legend

    def addLabels(self,drawing,title=None,xlabel=None,ylabel=None):
        from reportlab.graphics.charts.textlabels import Label
        if not title is None:
            Title = Label()
            Title.fontName = 'Helvetica-Bold'
            Title.fontSize = 7
            Title.x = drawing.width/2
            Title.y = drawing.height-25
            Title._text = title
            Title.maxWidth = 180
            Title.height = 20
            Title.textAnchor ='middle'
            drawing.add(Title)

        if not xlabel is None:
            XLabel = Label()
            XLabel.fontName = 'Helvetica'
            XLabel.fontSize = 7
            XLabel.x = drawing.width/2
            XLabel.y = 10
            XLabel.textAnchor ='middle'
            XLabel.maxWidth = 100
            XLabel.height = 20
            XLabel._text = xlabel
            drawing.add(XLabel)

        if not ylabel is None:
            YLabel = Label()
            YLabel.fontName = 'Helvetica'
            YLabel.fontSize = 7
            YLabel.x = 12
            YLabel.y = drawing.height/2
            YLabel.angle = 90
            YLabel.textAnchor ='middle'
            YLabel.maxWidth = 100
            YLabel.height = 20
            YLabel._text = ylabel
            drawing.add(YLabel)


    def colorRangeStatic(self,colorcount):
        """
        source: http://stackoverflow.com/questions/2328339/how-to-generate-n-different-colors-for-any-natural-number-n
        :param colorcount: Number of distinguishable colors to generate. Only pick first colorcount of 128, cycling
        :return: nothing
        """
        rgbhex = [
            "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
            "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
            "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
            "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
            "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
            "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
            "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
            "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

            "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
            "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
            "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
            "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
            "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
            "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
            "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
            "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"]

        self.colors = []
        for i in range(colorcount):
            ix = i%len(rgbhex)
            self.colors.append(colors.HexColor(rgbhex[ix]))


    def colorTest(self):
        from reportlab.graphics.shapes import Rect
        from reportlab.graphics.charts.textlabels import Label
        self.colorRangeStatic(130)

        dim=25
        width = self.PAGE_WIDTH-150
        inrow = 8 #int(width/dim)
        height = int(len(self.colors)/inrow)
        if len(self.colors)%inrow > 0:
            height +=1
        height *= dim
        drawing = Drawing(width, height)
        for i,col in enumerate(self.colors):
            x = (i%inrow)*dim
            y = int(i/inrow)*dim
            rec = Rect(x=x,y=y,width=dim,height=dim)
            rec.fillColor = col
            drawing.add(rec)
            lab = Label()
            lab.x=x+dim/2
            lab.y=y+dim/2
            lab.setText('%d'%i)
            drawing.add(lab)
        return drawing

    def stackedGraphs(self,lxydata,tmin=None,tmax=None): # list of (label, [x,yarr]) with x in seconds

        lcheight = 125
        lcwidth = self.PAGE_WIDTH-280
        graphsinchart = self.maxnumgraphsinchart
        stacks = (len(lxydata)/graphsinchart)
        if len(lxydata)%graphsinchart>0:
            stacks +=1
        drawing = Drawing(lcwidth+130, stacks*lcheight+75)
        ix = 0
        for ix in range(0,len(lxydata),graphsinchart):
            subset=lxydata[ix:min(ix+graphsinchart,len(lxydata))]
            maxX = subset[0][1][0][0]
            maxY = 1.
            minX = maxX
            minY = -1
            for label,xyarr in subset:
                xarr = [x for (x,y) in xyarr]
                yarr = [y for (x,y) in xyarr]
                maxY = max(maxY, max(yarr))
                minY = min(minY, min(yarr))
                maxX = max(maxX, max(xarr))
                minX = min(minX, min(xarr))
            if not tmin is None:
                minX = tmin
            if not tmax is None:
                maxX = tmax
            subset.append(('CritMax',[[minX,1.],[maxX,1.]]))
            subset.append(('CritMin',[[minX,-1.],[maxX,-1.]]))
            lc = LinePlot()
            lc.x = 30
            lc.y = 50+lcheight*ix/graphsinchart
            lc.height = (lcheight-10)
            lc.width = lcwidth
            lc.data = [xy for (l,xy) in subset]
            datalabels = [ l for (l,xy) in subset]

            lc.xValueAxis.valueMin = minX
            lc.xValueAxis.valueMax = maxX
            lc.xValueAxis.valueStep = int((maxX-minX)/5)
            if ix == 0:
                lc.xValueAxis.labelTextFormat = dateformatter # Use the formatter
            else:
                lc.xValueAxis.labelTextFormat = emptyformatter # Use the formatter
            lc.xValueAxis.visibleGrid           = True
            lc.xValueAxis.gridStrokeWidth = .5
            lc.xValueAxis.gridStrokeColor = colors.darkgrey
            lc.xValueAxis.subTickNum = 1
            lc.xValueAxis.subGridStrokeWidth = .25
            lc.xValueAxis.subGridStrokeColor = colors.lightgrey
            lc.xValueAxis.visibleSubGrid = True

            lc.yValueAxis.valueMin = minY
            lc.yValueAxis.valueMax = maxY
            lc.yValueAxis.valueStep = .25
            if (maxY-minY)>2:
                lc.yValueAxis.valueStep = .5
            if (maxY-minY)>4:
                lc.yValueAxis.valueStep = 1.
            lc.yValueAxis.visibleGrid           = True
            lc.yValueAxis.gridStrokeWidth = .5
            lc.yValueAxis.gridStrokeColor = colors.darkgrey
            lc.yValueAxis.subTickNum = 1
            lc.yValueAxis.subGridStrokeWidth = .25
            lc.yValueAxis.subGridStrokeColor = colors.lightgrey
            lc.yValueAxis.visibleSubGrid = True

            for x in range(len(lc.data)):
                lc.lines[x].strokeColor = self.colors[x]
                lc.lines[x].name = datalabels[x]
                if len(lc.data[x]) ==1:
                    lc.lines[x].symbol = makeMarker('FilledCircle') # added to make filled circles.
                lc.lines[x].strokeWidth = 2
                if lc.lines[x].name == 'CritMin' or lc.lines[x].name == 'CritMax': # distinguish min/max
                    lc.lines[x].strokeColor = self.colors[0]
                    lc.lines[x].strokeDashArray = (5, 1)

            self.addLabels(drawing,ylabel="geschaalde waarde")
            drawing.add(lc)
            ix += 1

            drawing.add(self.addAutoLegend(lc,(lc.x,lc.y,lc.width,lc.height),len(lc.data),side='right'))
        return drawing

