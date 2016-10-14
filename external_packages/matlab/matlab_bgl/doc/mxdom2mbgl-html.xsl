<?xml version="1.0" encoding="utf-8"?>

<!--
This is an XSL stylesheet which converts mscript XML files into HTML.
Use the XSLT command to perform the conversion.

Copyright 1984-2004 The MathWorks, Inc. 
$Revision$  $Date$
-->

<!DOCTYPE xsl:stylesheet [ <!ENTITY nbsp "&#160;"> <!ENTITY reg "&#174;"> ]>
<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:mwsh="http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd">
  <xsl:output method="html" indent="yes"/>
  <xsl:strip-space elements="mwsh:code"/>

<xsl:template match="mscript">
<html>
  <head>
<xsl:comment>
This HTML is auto-generated from an M-file.
To make changes, update the M-file and republish this document.
      </xsl:comment>
    <xsl:variable name="title">
      <xsl:value-of select="//steptitle[@style='document']"/>
    </xsl:variable>
    <title>
      <xsl:if test = "$title != ''"><xsl:value-of select="$title"/></xsl:if>
      <xsl:if test = "$title = ''"><xsl:value-of select="m-file"/></xsl:if>
    </title>
    <meta name="generator">
      <xsl:attribute name="content">MATLAB <xsl:value-of select="version"/></xsl:attribute>
    </meta>
    <meta name="date">
      <xsl:attribute name="content"><xsl:value-of select="date"/></xsl:attribute>
    </meta>
    <meta name="m-file">
      <xsl:attribute name="content"><xsl:value-of select="m-file"/></xsl:attribute>
    </meta>
    <link rel="stylesheet" type="text/css" href="../site.css" />
    <style>

body {
  background: white;
  color: black;
}

p.footer {
  text-align: right;
  font-size: xx-small;
  font-weight: lighter;
  font-style: italic;
  color: gray;
}

pre.codeinput {
  margin-left: 20px;
  margin-top: 10px;
  margin-bottom: 10px;
  background-color: #bbbbbb;
  border: solid 1px;
  font-size: 10pt;
  width: 620px;
}

p
{
	margin: 10px;
}

hr
{
    color: #bbbbbb;
    height: 4;
}

.main
{
	border-left-style: solid;
	margin-left: 100px;	
	width: 650px;
}

.upwhitesq
{
    position: relative;
    left: -5px;
    top: -8px;
    background: white;  
}
.downwhitesq
{
    position: relative;
    left: 95px;
    bottom: 10px;
    background: white;  
}

img
{
	text-align: center;
}

span.keyword {color: #0000FF}
span.comment {color: #228B22}
span.string {color: #A020F0}
span.untermstring {color: #B20000}
span.syscmd {color: #B28C00}

pre.showbuttons {
  margin-left: 30px;
  border: solid black 2px;
  padding: 4px;
  background: #EBEFF3;
}

pre.codeoutput {
  margin-left: 20px;
  margin-top: 10px;
  margin-bottom: 10px;
  font-size: 10pt;
  width: 520px;
}

pre.error {
  color: red;
}

.intro {
  width: 650px;
}

    </style>
  </head>
  <body>

    <!-- Determine if the there should be an introduction section. -->
    <xsl:variable name="hasIntro" select="count(cell[@style = 'overview'])"/>

    <!-- If there is an introduction, display it. -->
    <xsl:if test = "$hasIntro">
      <h1><xsl:value-of select="cell[1]/steptitle"/></h1>
      <introduction><div class="intro">
        <xsl:apply-templates select="cell[1]/text"/>
      </div></introduction>
    </xsl:if>
    
    <xsl:variable name="body-cells" select="cell[not(@style = 'overview')]"/>

    <!-- Include contents if there are titles for any subsections. -->
    <xsl:if test="count(cell/steptitle[not(@style = 'document')])">
      <xsl:call-template name="contents">
        <xsl:with-param name="body-cells" select="$body-cells"/>
      </xsl:call-template>
    </xsl:if>

    <div class="main">
    
    <!-- Loop over each cell -->
    <xsl:for-each select="$body-cells">
        <!-- Title of cell -->
        <xsl:if test="steptitle">
          
          <xsl:if test="position()&gt;1">
          <hr /><div class="upwhitesq">&nbsp;</div>
          </xsl:if>
          
          <xsl:variable name="headinglevel">
            <xsl:choose>
              <xsl:when test="steptitle[@style = 'document']">h1</xsl:when>
              <xsl:otherwise>h2</xsl:otherwise>
            </xsl:choose>
          </xsl:variable>
          <xsl:element name="{$headinglevel}">
            <xsl:value-of select="steptitle"/>
            <xsl:if test="not(steptitle[@style = 'document'])">
              <a>
                <xsl:attribute name="name">
                  <xsl:value-of select="position()"/>
                </xsl:attribute>
              </a>
            </xsl:if>
          </xsl:element>
        </xsl:if>

        <!-- Contents of each cell -->
        <xsl:apply-templates select="text"/>
        <xsl:apply-templates select="mcode-xmlized"/>
        <xsl:apply-templates select="mcodeoutput"/>
        <xsl:apply-templates select="img"/>

        
    </xsl:for-each>

    <hr /> <div class="upwhitesq">&nbsp;</div>
    </div>
    <div class="downwhitesq">&nbsp;</div>
    <!-- OLD FOOTER CODE <p class="footer">
      <xsl:value-of select="copyright"/><br/>
      Published with MATLAB&reg; <xsl:value-of select="version"/><br/>
    </p> -->
    
    <xsl:apply-templates select="originalCode"/>

  </body>
</html>
</xsl:template>


<xsl:template name="contents">
  <xsl:param name="body-cells"/>
  <h2>Contents</h2>
  <div><ul>
    <xsl:for-each select="$body-cells">
      <xsl:if test="./steptitle">        
        <li><a><xsl:attribute name="href">#<xsl:value-of select="position()"/></xsl:attribute><xsl:value-of select="steptitle"/></a></li>
      </xsl:if>
    </xsl:for-each>
  </ul></div>
</xsl:template>


<!-- HTML Tags in text sections -->
<xsl:template match="p">
  <p><xsl:apply-templates/></p>
</xsl:template>
<xsl:template match="ul">
  <div><ul><xsl:apply-templates/></ul></div>
</xsl:template>
<xsl:template match="li">
  <li><xsl:apply-templates/></li>
</xsl:template>
<xsl:template match="pre">
  <xsl:choose>
    <xsl:when test="@class='error'">
      <pre class="error"><xsl:apply-templates/></pre>
    </xsl:when>
    <xsl:otherwise>
      <pre><xsl:apply-templates/></pre>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>
<xsl:template match="b">
  <b><xsl:apply-templates/></b>
</xsl:template>
<xsl:template match="tt">
  <tt><xsl:apply-templates/></tt>
</xsl:template>
<xsl:template match="a">
  <a>
    <xsl:attribute name="href"><xsl:value-of select="@href"/></xsl:attribute>
    <xsl:apply-templates/>
  </a>
</xsl:template>

<!-- Code input and output -->

<xsl:template match="mcode-xmlized">
  <pre class="codeinput"><xsl:apply-templates/><xsl:text><!-- g162495 -->
</xsl:text></pre>
</xsl:template>

<xsl:template match="mcodeoutput">
    <pre class="codeoutput"><xsl:apply-templates/></pre>
</xsl:template>


<!-- ALTERNATIVE Code input -->

<xsl:template match="XXXmcode-xmlized">
  <a>
    <xsl:attribute name="name">c<xsl:value-of select="../mcode-count"/></xsl:attribute>
  </a>
  <pre class="showbuttons">
    <a>
      <xsl:attribute name="href">#c<xsl:value-of select="(../mcode-count)+1"/></xsl:attribute>
      <img align="right" src="http://www-internal.mathworks.com/images/clf/orange_arrow.gif" border="0"/>
    </a>
    <a>
      <xsl:attribute name="href">matlab:<xsl:value-of select="../mcode-flat"/></xsl:attribute>
      <img align="right" style="margin-left:4px;" src="http://www-internal/images/clf/yellow_arrow.gif
" border="0"/>
    </a>
    <xsl:apply-templates/>
  </pre>
</xsl:template>

<!-- Figure and model snapshots -->

<xsl:template match="img">
  <img vspace="5" hspace="5">
    <xsl:attribute name="src"><xsl:value-of select="@src"/></xsl:attribute><xsl:text> </xsl:text>
  </img>
</xsl:template>

<!-- Stash original code in HTML for easy slurping later. -->

<xsl:template match="originalCode">
  <xsl:variable name="xcomment">
    <xsl:call-template name="globalReplace">
      <xsl:with-param name="outputString" select="."/>
      <xsl:with-param name="target" select="'--'"/>
      <xsl:with-param name="replacement" select="'REPLACE_WITH_DASH_DASH'"/>
    </xsl:call-template>
  </xsl:variable>
<xsl:comment>
##### SOURCE BEGIN #####
<xsl:value-of select="$xcomment"/>
##### SOURCE END #####
</xsl:comment>
</xsl:template>

<!-- Colors for syntax-highlighted input code -->

<xsl:template match="mwsh:code">
  <xsl:apply-templates/>
</xsl:template>
<xsl:template match="mwsh:keywords">
  <span class="keyword"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:strings">
  <span class="string"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:comments">
  <span class="comment"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:unterminated_strings">
  <span class="untermstring"><xsl:value-of select="."/></span>
</xsl:template>
<xsl:template match="mwsh:system_commands">
  <span class="syscmd"><xsl:value-of select="."/></span>
</xsl:template>


<!-- Footer information -->

<xsl:template match="copyright">
  <xsl:value-of select="."/>
</xsl:template>
<xsl:template match="revision">
  <xsl:value-of select="."/>
</xsl:template>

<!-- Search and replace  -->
<!-- From http://www.xml.com/lpt/a/2002/06/05/transforming.html -->

<xsl:template name="globalReplace">
  <xsl:param name="outputString"/>
  <xsl:param name="target"/>
  <xsl:param name="replacement"/>
  <xsl:choose>
    <xsl:when test="contains($outputString,$target)">
      <xsl:value-of select=
        "concat(substring-before($outputString,$target),$replacement)"/>
      <xsl:call-template name="globalReplace">
        <xsl:with-param name="outputString" 
          select="substring-after($outputString,$target)"/>
        <xsl:with-param name="target" select="$target"/>
        <xsl:with-param name="replacement" 
          select="$replacement"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$outputString"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>
</xsl:stylesheet>

