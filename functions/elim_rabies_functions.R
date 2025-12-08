###############################################################################################################
#
# HEADER INFORMATION
#
# Project
#    Eliminate Rabies Trial
#
# Purpose
#    Provide functions to make tables and figures for the Eliminate Rabies Trial analysis
#
# Developed from functions originally written by
#    Alex McConnachie, Paul Johnson & Orla Doolin at the Robertson Centre for Biostatistics
#


##########################################################################################
#
# functions required while creating data set
#


# convert variable to exclude impossible values

excl.hilo<-
  function(x,Min=-Inf,Max=Inf)
  {
    x[(x<=Min)|(x>=Max)]<-NA
    x
  }

# function to create assumptions document and add new assumptions

assume<-
  function(assume.text,assume.num,assume.out,last=FALSE)
  {
    if(assume.num==1)
    {
      cat(
        "<html xmlns='http://www.w3.org/TR/REC-html40'>",
        '<head>',
        '</head>',
        '<body>',
        paste('<h1>',study.title,'</h1>',sep=''),
        paste('<h2>Assumptions Document</h2>',sep=''),
        '<br/>',
        paste('<p><em>Written by:</em><br/>',author,'</p>',sep=''),
        paste('<p><em>Analysis program last run on:</em><br/>',date(),' by ',run.by,'</p>',sep=''),
        paste('<p><em>Latest version of analysis program:</em><br/>',prog.file,'</p>',sep=''),
        paste('<p><em>Latest version of analysis report:</em><br/>',file.out,'</p>',sep=''),
        "<p><br/></p>",
        sep='\n',file=assume.out,append=FALSE)
    }
    cat(
      paste('<p><em>Assumption ',assume.num,'.</em><br/>',assume.text,'</p>',sep=''),sep='\n',file=assume.out,append=TRUE)
    if(last)
      cat("</body>","</html>",file=assume.out,sep="\n",append=TRUE)
  }


# function to show random sample of ordered rows of data frame, for inspecting data

sample.rows<-
  function(x,noof.rows=20,char.lim=NULL,dots=FALSE)
  {
    if(is.null(dim(x)))              # if x does not have >1 row and >1 column
    {
      smp<-x[sort(sample(1:length(x),noof.rows))]
      if(is.numeric(char.lim) & is.character(smp))
      {
        len.smp<-nchar(smp)
        smp<-substr(smp,1,char.lim)
        if(dots) smp[nchar(smp)<len.smp]<-paste(smp[nchar(smp)<len.smp],'...',sep='')
      }
    } else
    {
      smp<-x[sort(sample(1:nrow(x),noof.rows)),]
      if(is.numeric(char.lim))
        smp<-apply(smp,2,function(z) if(is.character(z)) substr(z,1,char.lim))
    }
    smp
  }


# function to split numberic vector around median in as even a way as possible

split5050<-
  function(x,less=TRUE)
  {
    med<-median(x,na.rm=TRUE)
    lt<-x<med
    le<-x<=med
    out<-if(abs(diff(table(lt))) < abs(diff(table(le)))) lt else le
    if(less) return(out)
    if(!less) return(!out)
  }

# function to set CRAN mirror, download (if required) and load package

install.load<-
  function(package) {
    if(!(package %in% installed.packages())) {
      mirrors<-getCRANmirrors()
      local.mirror<-mirrors$URL[mirrors$Country=="UK"][1]
      local({r <- getOption("repos"); r["CRAN"] <- local.mirror; options(repos=r)})
      install.packages(package)
    }
    eval(substitute(library(x),list(x=package)))
    return(paste("Package",package,"installed and loaded"))
  }

##########################################################################################
#
# Generic functions
#

# alex's SplusTables functions

# functions to detect if object is a table element
# allows vectors of sharacter or numeric as a Cell

is.Cell<-
  function(x)
  {
    if(is.vector(x)&(is.character(x)|is.numeric(x)))
      x<-list(Text=x)
    
    if(!is.list(x)|is.null(names(x)))
      return(F)
    
    !any(is.na(match(names(x),c("Text","Style","Other"))))
  }

is.Row<-
  function(x)
  {
    if(!is.list(x))
      return(F)
    
    all(unlist(lapply(x,is.Cell)))
  }

is.Table<-
  function(x)
  {
    if(!is.list(x))
      return(F)
    
    all(unlist(lapply(x,is.Row)))
  }

is.Element<-
  function(x)
    is.Cell(x)|is.Row(x)|is.Table(x)

# function to convert valid table element to table format

as.Table<-
  function(x)
  {
    if(!is.Element(x))
      stop("x is not a valid table element")
    
    if(is.vector(x)&(is.character(x)|is.numeric(x)))
      x<-list(list(list(Text=x)))
    
    if(is.Cell(x))
      x<-list(list(x))
    
    if(is.Row(x))
      x<-list(x)
    
    x
  }

# find number of rows or columns in table (assumes table rows are compatible)

NColsTable<-
  function(x)
    sum(unlist(lapply(as.Table(x)[[1]],function(x)max(1,x$Other$colspan,na.rm=T))))

NRowsTable<-
  function(x)
    length(as.Table(x))

# cbind function for lists
# if any element of list is a numeric or character vector, or a Cell or a Row, then convert to Table format
# if some Tables have 1 row and others have >1 row, then add rowspan=... to Tables with one row

cbindTable<-function(...)
{
  # convert arguments to a list
  args.list<-list(...)
  
  # check all arguments are table elements
  if(any(!unlist(lapply(args.list,is.Element))))
    stop("At least one argument is not a valid table element")
  
  # convert all arguments to Tables
  args.list<-lapply(args.list,as.Table)
  
  # check for compatibility
  # all must have 1 row, or the same number of rows
  args.rows<-unlist(lapply(args.list,NRowsTable))
  max.rows<-max(args.rows)
  compatible.rows<-(args.rows==1)|(args.rows==max.rows)
  args.id<-1:length(compatible.rows)
  if(!all(compatible.rows))
    stop(
      paste(
        "Incompatible number of rows in argument",
        if(sum(!compatible.rows)>1) "s",
        ": ",
        paste((args.id)[!compatible.rows],collapse=", ")))
  
  # add rowspan= to arguments with one column, if required
  if(max.rows>1)
    for(i in args.id[args.rows==1])
      args.list[[i]]<-
    lapply(
      args.list[[i]],
      function(xrow,max.rows)
      {
        xrow[[1]]$Other$rowspan<-max.rows
        xrow
      },
      max.rows=max.rows)
  
  # build the final table
  # run through i = 1:max.rows and pick out the i(th) row from each table and c() them together
  #   checks if a table has only 1 row, only the first is picked out;
  #   for later rows the "picking" function returns NULL, so that table does not contribute to the final c()'ed row
  
  lapply(
    as.list(1:max.rows),
    function(i,args.list)
      do.call(
        "c",
        lapply(
          args.list,
          function(x,i)
            if(length(x)>=i)
              x[[i]]
          else
            NULL,
          i=i)),
    args.list=args.list)
}


# rbind function for lists
# if any element of list is a numeric or character vector, or a Cell or a Row, then convert to Table format
# if some Tables have 1 column and others have >1 column, then add colspan=... to Tables with one column

rbindTable<-function(...)
{
  # convert arguments to a list
  args.list<-list(...)
  
  # check all arguments are table elements
  if(any(!unlist(lapply(args.list,is.Element))))
    stop("At least one argument is not a valid table element")
  
  # convert all arguments to Tables
  args.list<-lapply(args.list,as.Table)
  
  # check for compatibility
  # all must have 1 column, or the same number of columns
  args.cols<-unlist(lapply(args.list,NColsTable))
  max.cols<-max(args.cols)
  compatible.cols<-(args.cols==1)|(args.cols==max.cols)
  args.id<-1:length(compatible.cols)
  if(!all(compatible.cols))
    stop(
      paste(
        "Incompatible number of columns in argument",
        if(sum(!compatible.cols)>1) "s",
        ": ",
        paste((args.id)[!compatible.cols],collapse=", ")))
  
  # add colspan= to arguments with one column, if required
  if((max.cols>1)&any(args.cols==1))
    args.list<-
    lapply(
      args.list,
      function(xtable,max.cols)
      {
        if(NColsTable(xtable)==1)
          xtable<-
            lapply(
              xtable,
              function(xrow,max.cols)
              {
                xrow[[1]]$Other$colspan<-max.cols
                xrow
              },
              max.cols=max.cols)
        xtable
      },
      max.cols=max.cols)
  
  # build the final table
  # just c() them all together
  
  do.call("c",args.list)
}



# function to convert style info to text

StyleWrite<-
  function(Style=NULL,Other=NULL,inCell=T)
  {
    Other<-as.list(unlist(Other)[unique(names(Other))])
    Style<-as.list(unlist(Style)[unique(names(Style))])
    
    paste(
      paste(paste(names(Other),unlist(Other),sep="="),collapse=" "),
      if(!is.null(Style))
        paste(
          ifelse(inCell,"style='","{"),
          paste(paste(names(Style),unlist(Style),sep=":"),collapse=";"),
          ifelse(inCell,"'","}"),
          sep=""))
  }


# function to convert cell contents into text for writing to external file
# will use as do.call("CellWrite",CellData)

CellWrite<-
  function(Text="&nbsp",Style=NULL,Other=NULL)
  {
    Text<-
      ifelse(
        (Text=="")|((regexpr(" +",Text)==1)&(attr(regexpr(" +",Text),"match.length")==nchar(Text))),
        "&nbsp;",
        Text)
    paste(
      "<td ",StyleWrite(Style,Other),"><p style='page-break-after:avoid'>",
      #					"<td ",StyleWrite(Style,Other),"><p>",
      paste(Text,collapse="<br/>"),
      "</p></td>",
      sep="")
  }


# function to write a row to an external file

RowWrite<-
  function(x,AllCells=NULL,filename)
  {
    if(!is.Row(x))
      stop("x is not a row of a table")
    
    x<-
      lapply(
        x,
        function(y,AllCells)
        {
          y$Other<-c(y$Other,AllCells$Other)
          y$Style<-c(y$Style,AllCells$Style)
          y
        },AllCells=AllCells)
    
    cat(
      "<tr>",
      unlist(lapply(x,function(x)do.call("CellWrite",x))),
      "</tr>",
      file=filename,sep="\n",append=T)
  }


# export html table
# input
#   TableHeader - list to define top rows of table - title & column headings
#   TableRows - list to define data rows of table
#   AllCells - list of style & layout components common to all cells (unless superceded by individual cell definitions)
#   filename - output file
#
#   Each element of TableHeader and TableRows is a list defining one row of the table
#     Each row of table specified as a list, with each element defining a cell of the table
#
#   Each cell of table defined by:
#     list with elements
#       Text = cell contents - character vector, will have new lines inserted
#       Other = list of style & layout components for that cell e.g. alignment

ExportTable.HTML<-
  function(TableHeader,TableRows,AllCells=NULL,filename,first=F,last=F,fullWidth=T)
  {
    
    # if first table, then start document
    
    if(first)
      cat(
        "<html xmlns='http://www.w3.org/TR/REC-html40'>",
        "<head>",
        "</head>",
        "<body>",
        '<font face="Times New Roman">',
        sep="\n",file=filename,append=F)
    
    # if not first table, then add page break
    
    #			if(!first)
    #				cat(
    #					"<br style='page-break-before:always'>",
    #					sep="\n",file=filename,append=T)
    
    # start table
    
    cat(
      paste("<table",{if(fullWidth) " width=100%"},">",sep=""),
      file=filename,sep="\n",append=T)
    
    # write table header
    
    cat("<thead>",file=filename,sep="\n",append=T)
    
    lapply(TableHeader,RowWrite,AllCells=AllCells,filename=filename)
    
    cat("</thead>",file=filename,sep="\n",append=T)
    
    # write table data
    
    lapply(TableRows,RowWrite,AllCells=AllCells,filename=filename)
    
    # end table
    cat("</table>","<br/>",file=filename,sep="\n",append=T)
    
    # if last table, end html code
    if(last)
      cat("</body>","</html>",file=filename,sep="\n",append=T)
    
    invisible(NULL)
    
  }




# settings to allow ViewTable.HTML (below) to store html table in temporary directory, then delete it on exit

# set temporary directory

path.levels<- unlist(strsplit(tempfile(),"[\\]"))

.tmpdir<-
  paste(c(path.levels[1:max(1,(length(path.levels)-(2)))],'tmpStudyData'),collapse=ifelse(Sys.info()["sysname"]=="Windows","\\","/"))

if(!file.exists(.tmpdir)) dir.create(.tmpdir,recursive=TRUE)

# function to delete temporary directory on session end

.Last<-
  function(x)
    if(file.exists(.tmpdir)) unlink(.tmpdir,recursive=TRUE)


# function to allow viewing of html output in default browser

ViewTable.HTML<-
  function(xtable,legend.text="Legend text here.")
  {
    
    # create filename for table
    
    if(.Platform$OS.type=="windows") filename<-paste(tempfile(tmpdir=.tmpdir),'.html',sep='')
    if(.Platform$OS.type=="unix") filename<-"~/temp.table.html"
    print(filename)
    
    # use ExportTable.HTML() to create the table
    
    ExportTable.HTML(
      TableHeader=
        rbindTable(
          list(
            Text=paste("<b>Table X.</b> ",legend.text,sep=""),
            Other=list(align="left"),
            Style=list("border-bottom"="solid windowtext 1pt")),
          xtable[[1]]),
      TableRows=xtable[[2]],
      AllCells=
        list(
          Style=
            list(
              "font-size"="10.0pt",
              "border-collapse"="collapse","padding"=paste(rep("0.25cm",4),collapse=" "),
              "border-bottom"="solid windowtext 1pt",
              "border-top"="none","border-left"="none","border-right"="none"),
          Other=list(align="center",valign="middle")),
      filename=filename,
      first=TRUE,last=TRUE,fullWidth=FALSE)
    
    # open temp table file in internet explorer
    
    if(Sys.info()["sysname"]=="Windows") shell.exec(filename)
    if(substr(.Platform$pkgType,1,3)=="mac") system(paste("open /Applications/Firefox.app",filename))
  }


# function to output html table using alex's SplusTables functions

output.html.table<-
  function(
    xtable,legend.text,table.number=tabnum,file=file.out,first=FALSE,
    last=FALSE,fullWidth=FALSE,font.size=10,font.family="Times New Roman",landscape=FALSE,
    new.page=FALSE,new.page.after=FALSE,extra.break.after=TRUE)
  {
    legend.text<-
      paste('<b>Table ',table.number,'.&nbsp</b> ',legend.text,sep='')
    if(landscape)
      cat(
        "</div>",
        "<br clear=all style='page-break-before:always;mso-break-type:section-break'>",
        "<div class=Section2>",
        sep="\n",file=file,append=T)
    if(new.page & !landscape)
      cat(
        "<br style='page-break-after:always'>",
        sep="\n",file=file,append=T)
    ExportTable.HTML(
      TableHeader=
        rbindTable(
          list(
            Text=legend.text,
            Other=list(align="left"),
            Style=list("border-top"="none")),
          xtable$tab.heads),
      TableRows=xtable$tab.rows,
      AllCells=
        list(
          Style=
            list(
              "font-family"=font.family,
              "font-size"=paste(font.size,".0pt",sep=''),
              "border-collapse"="collapse","padding"=paste(rep("0.10cm",4),collapse=" "),
              "border-bottom"="solid windowtext 1pt",
              "border-top"="none","border-left"="none","border-right"="none"),
          Other=list(align="center",valign="middle")),
      filename=file,
      first=first,last=last,fullWidth=fullWidth)
    if(landscape)
      cat(
        "</div>",
        "<br clear=all style='page-break-before:always;mso-break-type:section-break'>",
        "<div class=Section1>",
        sep="\n",file=file,append=T)
    if(new.page.after & !landscape)
      cat(
        "<br style='page-break-after:always'>",
        sep="\n",file=file,append=T)
    if(extra.break.after)
      cat(
        "<br/>",
        sep="\n",file=file,append=T)
  }


# function to turn matrix into html-formatted table row

cbindMatrix<-
  function(x,use.rownames=FALSE,use.colnames=FALSE,Other=NULL,Style=NULL)
  {
    if(is.null(dim(x))) x<-matrix(x,ncol=length(x))
    if(use.colnames) x<-rbind(dimnames(x)[[2]],x)
    x.list<-
      lapply(as.list(1:ncol(x)),
             function(i,x,Other,Style) 
               cbindTable(list(Text=x[,i],Other=Other,Style=Style)),x=x,Other=Other,Style=Style)	
    out<-do.call("cbindTable",x.list)
    if(use.rownames) 
      out<-cbindTable(list(Text=dimnames(x)[[1]],Other=list(align="left")),out)
    out
  }



# function to set page orientation

landscape<-
  function(to.landscape,file)
  {
    section<-to.landscape+1
    cat(
      "</div>",
      "<br clear=all style='page-break-before:always;mso-break-type:section-break'>",
      paste("<div class=Section",section,">",sep=""),
      sep="\n",file=file,append=T)      
  }


# formatting

my.format<-
  function(x,ndp=0,na.string="")
  {
    out<-
      format(round(x,ndp),ns=ndp,scientific=FALSE,just='none')
    out[grep("NA",out)]<-na.string
    out
  }

# format p-values

p.format<-
  function(p,stars=FALSE,highlight=FALSE,p.eq=FALSE,alpha=0.05)
  {
    if(stars)
      p.formatted<-ifelse((p>alpha) | is.na(p),"",paste(rep("*",min(3,ceiling(-log10(p)-1),na.rm=TRUE)),collapse=""))
    else
    {
      p.formatted<-ifelse(as.numeric(my.format(p,3))<0.001,"<0.001",my.format(p,3))
      if(highlight) p.formatted<-ifelse(p<alpha,paste(p.formatted,"*",sep=""),p.formatted)
    }
    if(p.eq) p.formatted<-gsub("=<","<",paste("p=",p.formatted,sep=""))
    p.formatted[is.na(p)]<-"-"
    p.formatted
  }

# function to capitalise initial letters of words.
#   this function is adapted from the examples in the R Help for toupper()

capwords<-
  function(s,strict=FALSE,first.only=FALSE)
  {
    cap<-
      function(s)
        paste(
          toupper(substring(s,1,1)),
          {s <- substring(s,2); if(strict) tolower(s) else s},
          sep = "", collapse = " " )
    if(first.only)
      cap(s)
    else
      sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }


# make neater tables by deleting space after '('. also tidies up html output.

no.space.bracket<-
  function(x,html=TRUE)
  {
    x<-trim(gsub("  "," ",x))
    x<-gsub("  "," ",x)
    x<-gsub("NaN","-",x)
    x<-gsub('[\\(]  ','(',x)
    x<-gsub('[\\(] ','(',x)
    x<-gsub('[\\(][\\)]','(-)',x)
    x<-gsub('`','',x)
    if(html)
    {
      x<-gsub(' ','&nbsp;',x)
      x<-gsub('[\\+]','&#43;',x)
      x<-gsub(',','&#44;',x)
      x<-gsub('-','&#45;',x)
      x<-gsub('[\\.]','&#46;',x)
      x<-gsub(':','&#58;',x)
      x<-gsub('<','&lt;',x)
      x<-gsub('>','&gt;',x)
    }
    x
  }

# mean (SD)

mean.sd<-
  function(x)
  {
    mean<-mean(x,na.rm=TRUE)
    sd<-sd(x,na.rm=TRUE)
    paste(
      ifelse((is.na(mean) | is.nan(mean)),'-',my.format(mean,2)),' (',
      ifelse((is.na(sd) | is.nan(sd)),'-',my.format(sd,2)),')',sep='')
  }




###################################################################################################
#
#  Functions to create tables output
#


# function to output n(%)

n.perc<-
  function(x,names=FALSE,html=TRUE,ndp=0)
  {
    xtab<-table(x)
    xtab.names<-names(xtab)
    percentage<-100*xtab/sum(xtab)
    output.tab<-
      paste(
        xtab,' (',
        ifelse(
          (is.na(percentage) | is.nan(percentage)),
          '-)',
          paste(my.format(percentage,ndp),'%)',sep='')),
        sep='')
    if(html) out<-no.space.bracket(output.tab)
    names(out)<-xtab.names
    if(names) out<-cbind(no.space.bracket(xtab.names),out)
    out
  }



# function to output n(%) from table margins

n.perc.margin<-
  function(r.factor,c.factor,margin)
  {
    margin.tab<-margin.table(table(r.factor,c.factor),margin=margin)
    output<-
      paste(
        margin.tab,' (',
        my.format(100*margin.tab/sum(margin.tab),1),'%)',sep='')
    no.space.bracket(output)
  }


# function to output prevalence by cohort as n(%)

n.perc.prev<-
  function(cohort,outcome,test=TRUE)
  {
    xtab<-table(cohort,outcome)
    percentage<-100*xtab[,ncol(xtab)]/apply(xtab,1,sum)
    output<-
      c(
        paste(
          xtab[,ncol(xtab)],
          ' (',
          ifelse((is.na(percentage) | is.nan(percentage)),'-)',paste(my.format(percentage,1),'%)',sep='')),sep=''),
        {
          if(test)
            if(min(margin.table(xtab,margin=2))>0)
              p.format(my.fisher.test(xtab))
          else
            '-'
          else
            NULL
        })
    no.space.bracket(output)
  }

# function to put line breaks into long strings (useful for breaking up column headings

line.breaks<-
  function(string,max.char=20,break.string="\n")
  {
    break.str<-c("-")
    for(b.s in break.str)
      string<-gsub(b.s,paste(b.s," ",sep=""),string)
    words<-strsplit(string," ")
    word.lengths<-lapply(words,function(ww) nchar(ww)+1)
    cum.word.lengths<-lapply(word.lengths,cumsum)
    out<-
      unlist(
        lapply(
          1:length(words),
          function(w)
          {
            x<-cum.word.lengths[[w]]
            y<-x/max.char-floor(x/max.char)
            z<-y[-1]-y[-length(y)]<0
            line.no<-c(1,cumsum(z)+1)
            broken.lines<-
              lapply(
                1:max(line.no),
                function(l.n)
                  paste(words[[w]][line.no==l.n],collapse=" "))
            paste(unlist(broken.lines),collapse=break.string)
          }))
    breaks.loc<-
      gregexpr("\n",out)[[1]]
    attributes(out)<-
      list(lines=ifelse(length(breaks.loc)==1 && breaks.loc==-1,1,length(breaks.loc)+1))
    for(b.s in break.str)
      out<-gsub(paste(b.s," ",sep=""),b.s,out)
    gsub("</\n","</",gsub("/\n>","/>",out))
  }



# function to return ndp given a number and a desired number of sighnificant figures

ndp.count<-
  function(x,nsf=3)
  {
    if(x==0) return(nsf)
    max(floor(log10(1/abs(x)))+nsf,0,na.rm=TRUE)
  }

# function to summarise continuous variables and factors

my.summary<-
  function(x,text.format=TRUE,force.na.row=FALSE,ndp=2,auto.ndp=FALSE,miss.level=FALSE,sum.index=1:4,all.cats=TRUE) {
    x.na<-sum(is.na(x))
    x.na.n.perc<-n.perc(factor(is.na(x),c(T,F)))[1]
    names(x.na.n.perc)<-NULL
    if(is.numeric(x)) {
      if(length(na.omit(x))>1 && auto.ndp & sd(x,na.rm=TRUE)!=0) {
        ndp<-ndp.count(sd(x,na.rm=TRUE),ndp)
      }
      x<-na.omit(x)
      qx<-quantile(x)
      smry<-list(N=length(x),Mean=mean(x),'Std Dev'=sd(x),Median=median(x),Min=min(x),Max=max(x),Q1=qx[2],Q3=qx[4])
      if(text.format)
        smry<-
        list(
          'N<sub>OBS</sub> (N<sub>MISSING</sub>)'=paste(smry$N," (",x.na,")",sep=""),
          'Mean&nbsp;(SD)'=paste(my.format(smry$Mean,ndp)," (",my.format(smry$'Std Dev',ndp),")",sep=""),
          'Median&nbsp;(IQR)'=paste(my.format(smry$Median,ndp)," (",my.format(smry$Q1,max(ndp,0)),
                                    ", ",my.format(smry$Q3,max(ndp,0)),")",sep=""),
          '[Range]'=paste("[",my.format(smry$Min,max(ndp,0)),", ",my.format(smry$Max,max(ndp,0)),"]",sep=""))
      no.space.bracket(unlist(smry))[sum.index]
    } else {
      if(!is.factor(x)) x<-factor(x,levels(x))
      if(miss.level & any(is.na(x))) x<-factor(x,exclude=NULL,levels=c(levels(x),NA),labels=c(levels(x),"Missing"))
      levels(x)<-gsub("<","&lt;",levels(x))
      levels(x)<-gsub(">","&gt;",levels(x))
      out <- no.space.bracket(c("N<sub>OBS</sub> (N<sub>MISSING</sub>)"=paste(length(na.omit(x))," (",x.na,")",sep=""),n.perc(x)))[all.cats]
      if(!1 %in% sum.index) out <- out[-1]
      out
    }
  }




# do fisher exact on proportions or kruskal-wallis test on location of continuous variable

summary.tab.np.test<-
  function(x.data,test.cohort,simulate.p.value=TRUE)
  {
    if(length(unique(na.omit(test.cohort[!is.na(x.data)])))==1) return(NA)
    if(length(unique(na.omit(x.data[!is.na(test.cohort)])))==1) return(NA)
    if(is.factor(x.data)) return(fisher.test(table(x.data,test.cohort),simulate.p.value=simulate.p.value,B=10000)$p.value)
    if(!is.factor(x.data)) return(kruskal.test(as.numeric(x.data),test.cohort,correct=TRUE)$p.value)
  }

# example: 
# summary.tab.np.test(c(1,2,3,4,5),factor(c(1,1,1,0,0)))


# function to do summary stats of change in continuous variable by cohort

summary.tab<-
  function(x.data,cohort,cat.name=NULL,cohort.names=NULL,total=TRUE,n.cohort=FALSE,
           force.na.row=TRUE,n.cohort.na=TRUE,stat.head="",ndp=2,auto.ndp=FALSE,sum.index=1:4,all.cats=TRUE,test.list=NULL,
           p.below=FALSE)
  {
    if(is.numeric(x.data) && !is.factor(x.data) && length(unique(na.omit(x.data)))>1 && auto.ndp)
      ndp<-ndp.count(sd(x.data,na.rm=TRUE),ndp)
    if(is.null(cat.name)) cat.name<-label(cohort)
    if(n.cohort) levels(cohort)<-paste(levels(cohort),"<br/>(N=",table(cohort),")",sep="")
    if(is.null(cohort.names))
      cohort.names<-as.list(levels(cohort))
    cohort.na<-sum(is.na(cohort))
    if(n.cohort.na)
      cat.name<-paste(cat.name," (N<sub>MISSING</sub> = ",cohort.na,")",sep="")
    if(sum(is.na(cbind(x.data,cohort)))>0) force.na.row<-TRUE
    if(length(x.data)==0) # protect against difference between R 2.6.0 and later versions when splitting empty data frame
    {
      split.data<-lapply(levels(cohort),function(l) x.data)
      names(split.data)<-levels(cohort)
    } else
      split.data<-split(x.data,cohort)
    sum.stats<-
      lapply(split.data,my.summary,text.format=TRUE,force.na.row=force.na.row,ndp=ndp,sum.index=sum.index,all.cats=all.cats)
    tab.heads<-
      rbindTable(cat.name,do.call("cbindTable",cohort.names))
    if(total)
    {
      tab.heads<-cbindTable(ifelse(n.cohort,paste("Total<br/>(n=",sum(!is.na(cohort)),")",sep=""),"Total"),tab.heads)
      sum.stats<-
        c(total=list(my.summary(x.data,text.format=TRUE,force.na.row=force.na.row,ndp=ndp,sum.index=sum.index,all.cats=all.cats)),sum.stats)
    }
    
    tab.heads<-cbindTable(list(Text=stat.head,Other=list(align="left")),tab.heads)
    tab.rows.body<-do.call("cbindTable",lapply(sum.stats,as.vector))
    tab.rows.labels<-list(Text=names(sum.stats[[1]]),Other=list(align="left"))
    
    if(!is.null(test.list) & p.below)
    {
      if(is.logical(test.list) && test.list) test.list<-list(cohort)
      for(i in 1:length(test.list))
      {
        tc<-test.list[[i]]
        p.val<-summary.tab.np.test(x.data,tc)
        tab.rows.body<-rbindTable(tab.rows.body,p.format(p.val))
        tab.rows.labels<-rbindTable(tab.rows.labels,list(Text="P",Other=list(align="left")))
      }
    }
    
    tab.rows<-cbindTable(tab.rows.labels,tab.rows.body)
    
    if(!is.null(test.list) & !p.below)
    {
      if(is.logical(test.list) && test.list) {
        test.list<-list(cohort)
        for(i in 1:length(test.list))
        {
          tc<-test.list[[i]]
          p.val<-summary.tab.np.test(x.data,tc)
          tab.heads<-cbindTable(tab.heads,paste("P",names(test.list)[i],sep="<br/>"))
          tab.rows<-cbindTable(tab.rows,p.format(p.val))
        }
      }
    }
    list(tab.heads=tab.heads,tab.rows=tab.rows)
  }



# function to rbind together summary.tab

demog.tab<-
  function(x.data,cohort,x.data.names=NULL,cat.name=NULL,total=TRUE,n.cohort=FALSE,n.cohort.na=FALSE,
           force.na.row=TRUE,allow.na.col=FALSE,stat.head="",footnote=FALSE,ndp=2,auto.ndp=TRUE,sum.index=1:4,
           all.cats=TRUE,test.list=NULL,p.below=FALSE,max.char=20)
  {
    
    if(allow.na.col & any(is.na(cohort)))
      cohort<-na.factor(cohort)
    
    if(is.null(x.data.names))
      x.data.names<-
        unlist(
          lapply(
            x.data,
            function(x)
              if(!is.null(attr(x,"subset")))
                paste(label(x)," (N&nbsp;=&nbsp;",sum(attr(x,"subset"),na.rm=TRUE),")",sep="") else
                  label(x)))
    
    if(is.null(cat.name))
      cat.name<-label(cohort)
    
    tab.main<-
      lapply(
        x.data,
        function(x,cht,cat.name,total,force.na.row,stat.head,ndp,auto.ndp,sum.index,all.cats,test.list,p.below)
        {
          if(!is.null(attr(x,"subset")))
          {
            cht<-cht[attr(x,"subset")]
            x<-x[attr(x,"subset")]
          }
          summary.tab(x.data=x,cat.name=cat.name,cohort=cht,total=total,force.na.row=force.na.row,
                      stat.head=stat.head,ndp=ndp,auto.ndp=auto.ndp,n.cohort=n.cohort,n.cohort.na=n.cohort.na,
                      sum.index=sum.index,all.cats=all.cats,test.list=test.list,p.below=p.below)
        },
        cht=cohort,cat.name=cat.name,total=total,force.na.row=force.na.row,
        stat.head=stat.head,ndp=ndp,auto.ndp=auto.ndp,sum.index=sum.index,all.cats=all.cats,test.list=test.list,
        p.below=p.below)
    
    
    tab.heads<-
      cbindTable(list(Text="",Other=list(align="left")),tab.main[[1]]$tab.heads)
    tab.rows<-        
      do.call(
        "rbindTable",
        lapply(
          1:length(tab.main),
          function(i,tab.main,x.data.names)
          {
            cbindTable(
              list(
                Text=line.breaks(x.data.names[i],break.string="<br/>",max.char=max.char),
                Other=list(valign="top",align="left")),
              tab.main[[i]]$tab.rows)
          },tab.main=tab.main,x.data.names=x.data.names))
    
    if(!allow.na.col & footnote)
      tab.rows<-
      rbindTable(tab.rows,
                 list(
                   Text=paste("Note: Missing data in the column variable was excluded",
                              "from the subgroup columns but included in the total column.",
                              "The subgroups columns will therefore not always sum to the total column."),
                   Other=list(align="left")))
    list(tab.heads=tab.heads,tab.rows=tab.rows)
  }

# function to cbind together output of demog.tab

demog.tab.wide<-
  function(x.data,cohort.list)
  {
    tab.list<-
      lapply(cohort.list,
             function(cohort) demog.tab(x.data=x.data,cohort=cohort,total=FALSE,test.list=FALSE,p.below=TRUE))
    row.head<-
      do.call("rbindTable",
              lapply(tab.list[[1]]$tab.rows,
                     function(x)
                     { 
                       if(x[[1]]$Text=="P") x[1] else x[1:2]
                     }))
    for(j in 1:length(row.head)) row.head[[j]][[1]]$Other$colspan<-1
    NRowsTable(row.head)
    row.body<-
      do.call("cbindTable",
              lapply(tab.list,
                     function(x)
                       do.call("rbindTable",
                               lapply(x$tab.rows,function(xx)cbindTable("",do.call("cbindTable",if(xx[[1]]$Text=="P") xx[2] else xx[-1:-2]))))))
    tab.heads.raw<-
      do.call("cbindTable",
              lapply(tab.list,
                     function(x)cbindTable("",rbindTable(x$tab.heads[[1]][-1:-2],x$tab.heads[[2]]))))
    tab.heads<-cbindTable(tab.list[[1]]$tab.heads[[1]][1:2],tab.heads.raw)
    NRowsTable(row.body)
    tab.rows<-cbindTable(row.head,row.body)
    list(tab.heads=tab.heads,tab.rows=tab.rows)
  }




# SUMMARY FUNCTION FOR LOGISTIC REGRESSION MODEL
lr.summ<-
  function(fit0,term=1,pval=T,ndp=1,trans=exp,p.only=F,sig.bold=0,p.break=F,fit1=NULL,...)
  {
    t.df<-fit0$df.residual
    t.thresh<-qt(0.975,t.df)
    coef0<-summary(fit0)$coef[term+1,,drop=F]
    
    if(!p.only & is.null(fit1))
    {
      p<-coef0[,4]
      out<-
        paste(
          my.format(trans(coef0[,1]),ndp),
          " (",my.format(trans(coef0[,1]-t.thresh*coef0[,2]),ndp),", ",my.format(trans(coef0[,1]+t.thresh*coef0[,2]),ndp),")",sep="")
      if(pval)
        if(p.break)
          out<-paste(out,"<br>",p.format(p,...),sep="")
      else
        out<-paste(out,", ",p.format(p,...),sep="")
      out
    }
    
    else
    {
      p<-if(is.null(fit1)) coef0[,4] else na.omit(anova(fit0,fit1,test="Chisq")[,"P(>|Chi|)"])
      out<-p.format(p,...)
    }
    
    out<-ifelse(p<sig.bold,paste("<b>",out,"</b>",sep=""),out)
    
    gsub("\\( ","\\(",as.vector(out))
    
  }


# SUMMARY FUNCTION FOR Firth's penalized-likelihood LOGISTIC REGRESSION MODEL 
# using the logistf package
# This function differs from lr.summ in not having an anova option, so there is no
# "fit1" argument
#
# for more details of Firth's penalized-likelihood logistic regression
# see:
# Statistical Methods in Medical Research 2008; 17: 621?634
# A comparative study of the bias corrected estimates in logistic regression
# Tapabrata Maiti & Vivek Pradhan
# www.ncbi.nlm.nih.gov/pubmed/18375454  
# 

logistf.summ<-
  function(fit0,term=1,pval=T,ndp=1,trans=exp,p.only=F,sig.bold=0,p.break=F,...)
  {
    coef0<-coef(fit0)[term+1]
    p<-fit0$prob[term+1]
    
    if(!p.only)
    {
      out<-
        paste(
          my.format(trans(coef0),ndp),
          " (",my.format(trans(fit0$ci.lower[term+1]),ndp),", ",my.format(trans(fit0$ci.upper[term+1]),ndp),")",sep="")
      if(pval)
        out<-if(p.break) paste(out,"<br>",p.format(p,...),sep="")	else paste(out,", ",p.format(p,...),sep="")
      out
    }
    
    else
    {
      out<-p.format(p,...)
    }
    
    out<-ifelse(p<sig.bold,paste("<b>",out,"</b>",sep=""),out)
    
    gsub("\\( ","\\(",as.vector(out))
    
  }


# SUMMARY FUNCTION FOR LOGISTIC REGRESSION MODEL WITH EXACT INFERENCE.
# AN IMPORTANT DIFFERENCE FROM lr.summ IS THAT THERE IS NO term ARGUMENT.
# IT'S ASSUMED THAT THE TERM OF INTEREST HAS BEEN SELECETED IN THE UPSTREAM
# FUNCTIONS (elrm OR elrm.fit) USING THE interest ARGUMENT.

elr.summ<-
  function(fit0,pval=T,ndp=1,trans=exp,p.only=F,sig.bold=0,p.break=F,...)
  {
    coef0<-fit0$coeffs
    
    if(!p.only)
    {
      p<-fit0$p.values
      out<-
        paste(
          my.format(trans(coef0),ndp),
          " (",my.format(trans(fit0$coeffs.ci$lower),ndp),", ",my.format(trans(fit0$coeffs.ci$upper),ndp),")",sep="")
      if(pval)
        out<-ifelse(p.break,paste(out,"<br>",p.format(p,...),sep=""),paste(out,", ",p.format(p,...),sep=""))
      out
    }
    
    else
    {
      p<-fit0$p.values
      out<-p.format(p,...)
    }
    
    out<-ifelse(p<sig.bold,paste("<b>",out,"</b>",sep=""),out)
    
    gsub("\\( ","\\(",as.vector(out))
    
  }


# SUMMARY FUNCTION FOR LOGISTIC REGRESSION MODEL  - p value output only
lr.summ.pval<-
  function(fit0,term=1,pval=T,ndp=1,trans=exp,p.only=F,sig.bold=0,p.break=F,fit1=NULL,...)
  {
    t.df<-fit0$df.residual
    t.thresh<-qt(0.975,t.df)
    coef0<-summary(fit0)$coef[term+1,,drop=F]
    
    if(!p.only & is.null(fit1))
    {
      p<-coef0[,4]
      out<-
        paste(
          my.format(trans(coef0[,1]),ndp),
          " (",my.format(trans(coef0[,1]-t.thresh*coef0[,2]),ndp),", ",my.format(trans(coef0[,1]+t.thresh*coef0[,2]),ndp),")",sep="")
      if(pval)
        if(p.break)
          out<-p.format(p,...)
      else
        out<-p.format(p,...)
      out
    }
    
    else
    {
      p<-if(is.null(fit1)) coef0[,4] else na.omit(anova(fit0,fit1,test="Chisq")[,"P(>|Chi|)"])
      out<-p.format(p,...)
    }
    
    out<-ifelse(p<sig.bold,paste("<b>",out,"</b>",sep=""),out)
    
    gsub("\\( ","\\(",as.vector(out))
    
  }




# function to show N (%) prevalence of the 2nd level of a two-level factor, or mean (SD) of a numeric variable,
# within subgroups of column factors, overall and by row factors
# 
# odds ratios with 95% CI and p-value are produced. by default these are are estimated via glm(...,family="binomial").
# these estimates are likely to be biased or undefined at with low cell counts. in this case, exact logistic regression 
# can be performed via elrm by supplying elrm.pars (see elrm.fit for details). for example:
# elrm.pars=list(round.ndp=0,burnin=5000,iter.used=100000,iter.ext=10000,max.updates=Inf,max.p.value.se=0,min.effective.iter=10000,r=2)
# will run until 10000 effective MCMC samples have been saved, or until muemory runs out.
# alternatively, Firth's penalized-likelihood logistic regression (www.ncbi.nlm.nih.gov/pubmed/18375454)
# can be used by setting firth.lrm=TRUE. This will over-ride elrm.pars.
# 

risk.factor.tab<-
  function(predictor,row.list,col.list,x.data,ndp=1,row.list.adjust=FALSE,subgroups=TRUE,firth.lrm=FALSE,elrm.pars=NULL)
  {
    if(
      (!is.numeric(x.data[,predictor]) & !is.factor(x.data[,predictor])) |
      (is.factor(x.data[,predictor]) && length(levels(x.data[,predictor]))!=2)) stop("The predictor must be a number or a factor with 2 levels")
    if(!is.null(dim(apply(sapply(x.data[,col.list],levels),1,unique)))) stop("Column factor must all have the same levels in the same order")
    tab.heads<-cbindTable("","",do.call("cbindTable",lapply(x.data[,col.list],function(zzz)c(line.breaks(label(zzz),break.string="<br/>")))))
    if(subgroups) tab.heads<-cbindTable("",tab.heads)
    tab.rows<-  
      do.call("rbindTable",lapply(0:ifelse(subgroups,length(row.list),0),
                                  function(i)
                                  {
                                    if(i==0)
                                    {
                                      z<-rep("",nrow(x.data))
                                      label(z)<-"All"
                                    } else 
                                      z<-x.data[,row.list[i]]
                                    x<-
                                      do.call("cbindTable",as.list(lapply(col.list,
                                                                          function(y)
                                                                          {
                                                                            tab.by.rowfactor<-
                                                                              by(x.data[,c(predictor,y,row.list)],z,
                                                                                 function(xy)
                                                                                 {
                                                                                   mod.form<-paste("~",predictor)
                                                                                   if(row.list.adjust & i==0) mod.form<-paste(mod.form,"+",paste(row.list,collapse="+"))
                                                                                   if(is.null(elrm.pars) | firth.lrm)
                                                                                   {
                                                                                     mod.form.full<-paste("I(",y,"=='",levels(x.data[,y])[2],"')",mod.form,sep="")
                                                                                     mod<-
                                                                                       if(firth.lrm) logistf(formula(mod.form.full),data=xy) else 
                                                                                         glm(mod.form.full,family="binomial",control=list(maxit=25),data=xy)
                                                                                   } else
                                                                                   {
                                                                                     mod<-
                                                                                       elrm.fit(outcome=y,outcome.level=levels(xy[,y])[2],formula.rhs=formula(mod.form),
                                                                                                interest=formula(paste("~",predictor)),source.data=xy,round.ndp=elrm.pars$round.ndp,burnin=elrm.pars$burnin,
                                                                                                iter.used=elrm.pars$iter.used,iter.ext=elrm.pars$iter.ext,max.updates=elrm.pars$max.updates,
                                                                                                max.p.value.se=elrm.pars$max.p.value.se,min.effective.iter=elrm.pars$min.effective.iter,
                                                                                                r=elrm.pars$r)
                                                                                   }
                                                                                   mod.summ<-
                                                                                     if(is.null(elrm.pars) | firth.lrm)
                                                                                     {
                                                                                       if(firth.lrm)
                                                                                       {
                                                                                         logistf.summ(fit0=mod,term=1,p.break=TRUE,p.eq=TRUE,p.only=FALSE,ndp=ndp)
                                                                                       } else
                                                                                       {
                                                                                         ifelse(mod$converged,lr.summ(fit0=mod,term=1,p.break=TRUE,p.eq=TRUE,p.only=FALSE,ndp=ndp),"-")
                                                                                       }
                                                                                     } else
                                                                                     {
                                                                                       elr.summ(fit0=mod,p.break=TRUE,p.eq=TRUE,p.only=FALSE,ndp=ndp)
                                                                                     }
                                                                                   
                                                                                   out<-
                                                                                     c(
                                                                                       by(xy[,predictor],xy[,y],
                                                                                          function(w)
                                                                                            paste(
                                                                                              no.space.bracket(
                                                                                                my.summary(w,all.cats=3,sum.index=2,ndp=ndp)),
                                                                                              "<br/>n=",sum(!is.na(w)),sep="")),
                                                                                       mod.summ)
                                                                                   
                                                                                   
                                                                                   out.tab<-cbindTable(out)
                                                                                   if(y==col.list[1])
                                                                                   {
                                                                                     names(out)[1]<-paste(names(out)[1],"<br/>",sep="")
                                                                                     out.tab<-
                                                                                       cbindTable(
                                                                                         cbindTable(
                                                                                           list(Text=c(names(out),"OR<br/>&nbsp;"),
                                                                                                Other=list(align="left"))),out.tab)
                                                                                   }
                                                                                   out.tab
                                                                                 })
                                                                            out<-do.call("rbindTable",as.list(tab.by.rowfactor))
                                                                            if(y==col.list[1])
                                                                            {
                                                                              out<-
                                                                                cbindTable(
                                                                                  do.call("rbindTable",
                                                                                          lapply(names(tab.by.rowfactor),function(nn)list(Text=nn,Other=list(align="left")))),
                                                                                  out)
                                                                            }
                                                                            out
                                                                          })))
                                    if(subgroups) cbindTable(list(Text=label(z),Other=list(align="left")),x) else x
                                  }))
    list(tab.heads=tab.heads,tab.rows=tab.rows)
  }






# function to produce table of factor analysis loadings
#
# added method="spearman" for rank correlations into a FA
#
# added switch to do factor analysis on polychoric and polyserial 
# correlations using hetcor {polycor} IF any of the variables 
# are *ordered* factors
# WARNING - this only works if the variables in the data frame 
# haven't been labelled. I suggest making a new data frame 
# of the numeric variables with no labels, before converting 
# to ordered factors (because the label attribute prevents
# hetcor from seeing the "levels" and "class" attributes.

loadings.tab<-
  function(
    data.list,n.factors,bold.rows=NULL,sorted=TRUE,gapped=TRUE,rotation="promax",start=NULL,threshold=0.4,
    table.out=TRUE,plot=FALSE,simu=50,boot=NULL,ndp=2,method="pearson")
  {
    if(any(sapply(data.list$data,is.ordered)))
    {
      require(polycor)
      print("Using polychoric correlations")
      cor.fun<-function(x)hetcor(x,std.err=FALSE,ML=FALSE,use="pairwise.complete",pd=TRUE)$correlations
    } else
    {
      cor.fun<-function(x)cor(x,use="pairwise.complete.obs",method=method)
    }
    data.names<-
      if(is.null(data.list$names)) {sapply(data.list$data,label)} else {data.list$names} 
    factanal.data<-data.list$data
    #          as.data.frame(lapply(data.list$data,function(x)if(is.factor(x))as.numeric(x)else x))
    factanal.data<-factanal.data[apply(!is.na(factanal.data),1,sum)>=2,]          
    na.pc<-paste(round(100*mean(is.na(factanal.data))),"%",sep="")
    n<-nrow(factanal.data)
    if(plot)neat.scree.plot(factanal.data,simu=simu,method=method)
    if(table.out)
    {
      
      fa.out<-
        factanal(covmat=cor.fun(factanal.data),factors=n.factors,rotation=rotation,start=start)
      load.tab<-fa.out$loadings[,]
      if(is.null(dim(load.tab))) load.tab<-cbind(Factor1=load.tab)
      if(threshold!=0)
        load.tab<-load.tab*(abs(load.tab)>=threshold)
      num.load.tab<-load.tab  # store numeric version and use later to order rows
      num.load.tab<-t(apply(num.load.tab,1,function(x){x[abs(x)<max(abs(x))]<-0;x}))
      if(n.factors==1) num.load.tab<-t(num.load.tab)
      load.tab<-as.matrix(my.format(load.tab,ndp=ndp))
      if(is.numeric(boot))
      {
        start.time<-Sys.time()
        boot.quantiles<-
          boot.loadings(factanal.data,n.factors=n.factors,rotation=rotation,R=boot,quant.list=list(0.025,0.975))
        end.time<-as.numeric(Sys.time())-as.numeric(start.time)
        print(
          paste(
            boot,'bootstrap replicates took',
            my.format(floor(end.time/60)),'min and',
            my.format(end.time-floor(end.time/60)*60),'s.'))
        if(boot.quantiles[[3]]<=0.1*boot)
        {
          for(i in 1:ncol(load.tab))
          {
            load.tab[,i]<-
              paste(
                load.tab[,i],
                ' (',
                my.format(boot.quantiles[[1]][,i],ndp=ndp),
                ', ',
                my.format(boot.quantiles[[2]][,i],ndp=ndp),')',
                sep='')
          } 
        }
        else 
        {
          warning('Confidence intervals not estimated because >10% of bootstrap replicates failed.',call.=FALSE)
          boot<-NULL
        }
      }
      if(threshold!=0)
        load.tab[regexpr(paste(c(0,'.',rep(0,ndp)),collapse=''),load.tab)>0]<-'.'
      load.tab<-cbind(data.names,no.space.bracket(load.tab))
      if(sorted)
      {
        for(i in max(2,ncol(load.tab)):2)
        {
          load.tab<-load.tab[order(abs(num.load.tab[,i-1]),decreasing=TRUE),]
          num.load.tab<-num.load.tab[order(abs(num.load.tab[,i-1]),decreasing=TRUE),]
          if(is.null(dim(num.load.tab))) num.load.tab<-cbind(num.load.tab)
        }
        if(gapped)  
        {
          load.tab<-
            do.call("rbind",
                    lapply(
                      1:nrow(load.tab),
                      function(i)
                      {
                        if(i %in% (1:nrow(num.load.tab))[apply(apply(num.load.tab==0,2,diff)==1,1,sum)==1])
                          rbind(load.tab[i,],"")
                        else load.tab[i,]
                      }))
        }
      }
      if(!is.null(bold.rows))
        load.tab[bold.rows,]<-gsub("<b></b>","",paste("<b>",load.tab[bold.rows,],"</b>",sep=""))
      factor.names<-gsub('Factor','Factor ',colnames(load.tab)[-1])
      tab.heads<-
        cbindTable(
          list(Text='Item',Other=list(align='left')),
          do.call("cbindTable",as.list(factor.names)))
      tab.rows<-
        do.call('cbindTable',
                lapply(
                  1:ncol(load.tab),
                  function(ii)list(Text=load.tab[,ii],Other=list(align=ifelse(ii==1,'left','center')))))
      notes<-
        paste(
          "Note: ",
          ifelse(n.factors>1,paste('Loadings were estimated using the <em>',rotation,'</em> rotation. ',sep=""),""),
          ifelse(threshold>0,paste('Only loadings exceeding &plusmn; ',threshold,' are shown.',sep=''),''),
          ifelse(
            is.numeric(boot),
            paste(' Confidence intervals were estimated from',boot-boot.quantiles[[3]],'bootstrap replicates.'),
            ''),
          sep='')
      tab.rows<-
        rbindTable(tab.rows,list(Text=notes,Other=list(align='left')))
      
      # make table of factor correlations
      
      if(!is.null(fa.out$rotmat)) 
      {
        tmat<-solve(fa.out$rotmat)
        fact.cor<-tmat %*% t(tmat)
        fact.cor<-cbind(factor.names,my.format(fact.cor,2))
        fact.cor.tab<-
          list(
            tab.heads=
              cbindTable(
                list(Text='',Other=list(align='left')),
                do.call("cbindTable",as.list(factor.names))),
            tab.rows=
              do.call('cbindTable',
                      lapply(
                        1:ncol(fact.cor),
                        function(j)list(Text=fact.cor[,j],Other=list(align=ifelse(j==1,'left','center'))))))                
      } else fact.cor.tab<-NULL
      list(tab.heads=tab.heads,tab.rows=tab.rows,fact.cor.tab=fact.cor.tab,load.tab=load.tab,fa.out=fa.out,n=n,na.pc=na.pc)
    }
  }


# function to produce bootstrapped stats from loadings

boot.loadings<-
  function(x.data.frame,n.factors,rotation="promax",R=100,quant.list=list(0.025,0.975))
  {
    
    # choose correlation function
    
    if(any(sapply(x.data.frame,is.ordered)))
    {
      require(polycor)
      print("Bootstrapping using polychoric correlations")
      cor.fun<-function(x)hetcor(x,std.err=FALSE,ML=FALSE,use="pairwise.complete",pd=TRUE)$correlations
    } else
    {
      cor.fun<-function(x)cor(x,use="pairwise.complete.obs",method=method)
    }
    
    # produce bootstrapped data.sets
    
    boot.data<-
      lapply(
        1:R,
        function(i)
          x.data.frame[sample(nrow(x.data.frame),replace=TRUE),])
    
    # factor analysis without bootstrap sets factor order for bootstapped data
    # first do the factor analysis then run function that gets name of variable
    # contributing most to each factor
    
    fa<-
      as.matrix(factanal(covmat=cor.fun(x.data.frame),factors=n.factors,rotation=rotation)$loadings[,])
    get.fact.names<-
      function(loadings)
      {
        unlist(lapply(
          1:ncol(loadings),
          function(c)
          {
            is.max.ld<-
              abs(loadings[,c])==apply(loadings,2,function(qq)max(abs(qq)))[c]
            is.max.ld<-is.max.ld[is.max.ld]
            names(is.max.ld)
          }))
      }
    
    fact.names<-
      get.fact.names(fa)
    
    # do factor analysis on bootstrapped data sets
    boot.load<-
      lapply(
        boot.data,
        function(bd)
          as.matrix(factanal(covmat=cor.fun(bd),factors=n.factors,rotation=rotation)$loadings[,]))
    
    # re-order columns of boot.load to match the fa loadings.
    # the factor names were defined above as the name of the variable with the highest loading.
    # now the highest loading for that variable is assumed to identify the same factor.
    # it is possible for a variable to be the highest in one factor but have a higher loading in
    # another factor where it isn't the highest. this will result in a warning, and if it happens frequently it probably
    # suggests that either n is too low or too many factors are being assumed.
    
    ordered.boot.load<-  
      do.call('cbind',lapply(
        boot.load,
        function(bl)
        {        
          col.order<-
            unlist(lapply(fact.names,function(f.n)order(abs(bl[f.n,]),decreasing=TRUE)[1]))
          if(length(unique(col.order))<length(col.order)){
            bl<-bl+NA
          }
          return(bl[,col.order])
        }))
    order.fail<-sum(is.na(ordered.boot.load))/(nrow(ordered.boot.load)*n.factors)
    if(order.fail>1){
      warning(
        paste(
          'Re-ordering of loadings columns failed in ',
          order.fail,
          ' bootstrap replicates (',
          my.format(100*order.fail/R,1),
          '%). High failure rate may indicate that you are assuming too many factors.',
          sep=''))}
    
    # produce desired statistic from replicated data sets
    
    output<-
      lapply(
        quant.list,
        function(quant){
          do.call('cbind',lapply(
            as.list(1:n.factors),
            function(j){
              apply(ordered.boot.load[,(1:R-1)*n.factors+j],1,function(nn)quantile(nn,probs=quant,na.rm=TRUE))
            }))})
    names(output)<-paste('q',1:length(output),sep='')
    c(output,order.fail)
  }




##use varimax rotation
# added method="spearman" for rank correlations into a FA
loadings.tab.varimax<-
  function(
    data.list,n.factors,bold.rows=NULL,sorted=TRUE,gapped=TRUE,rotation="varimax",start=NULL,threshold=0.4,
    table.out=TRUE,plot=FALSE,simu=50,boot=NULL,ndp=2, method='pearson')
  {
    data.names<-
      if(is.null(data.list$names)) {sapply(data.list$data,label)} else {data.list$names} 
    factanal.data<-
      as.data.frame(lapply(data.list$data,function(x)if(is.factor(x))as.numeric(x)else x))
    factanal.data<-factanal.data[apply(!is.na(factanal.data),1,sum)>=2,]          
    na.pc<-paste(round(100*mean(is.na(factanal.data))),"%",sep="")
    n<-nrow(factanal.data)
    if(plot)neat.scree.plot(factanal.data,simu=simu)
    if(table.out)
    {
      fa.out<-
        factanal(covmat=cor(factanal.data,use='pairwise.complete.obs',method=method),factors=n.factors,rotation=rotation,start=start)
      load.tab<-fa.out$loadings[,]
      if(threshold!=0)
        load.tab<-load.tab*(abs(load.tab)>=threshold)
      num.load.tab<-load.tab  # store numeric version and use later to order rows
      num.load.tab<-t(apply(num.load.tab,1,function(x){x[abs(x)<max(abs(x))]<-0;x}))          
      load.tab<-as.matrix(my.format(load.tab,ndp=ndp))
      if(is.numeric(boot))
      {
        start.time<-Sys.time()
        boot.quantiles<-
          boot.loadings(factanal.data,n.factors=n.factors,rotation=rotation,R=boot,quant.list=list(0.025,0.975))
        end.time<-as.numeric(Sys.time())-as.numeric(start.time)
        print(
          paste(
            boot,'bootstrap replicates took',
            my.format(floor(end.time/60)),'min and',
            my.format(end.time-floor(end.time/60)*60),'s.'))
        if(boot.quantiles[[3]]<=0.1*boot)
        {
          for(i in 1:ncol(load.tab))
          {
            load.tab[,i]<-
              paste(
                load.tab[,i],
                ' (',
                my.format(boot.quantiles[[1]][,i],ndp=ndp),
                ', ',
                my.format(boot.quantiles[[2]][,i],ndp=ndp),')',
                sep='')
          } 
        }
        else 
        {
          warning('Confidence intervals not estimated because >10% of bootstrap replicates failed.',call.=FALSE)
          boot<-NULL
        }
      }
      if(threshold!=0)
        load.tab[regexpr(paste(c(0,'.',rep(0,ndp)),collapse=''),load.tab)>0]<-'.'
      load.tab<-cbind(data.names,no.space.bracket(load.tab))
      if(sorted)
      {
        for(i in ncol(load.tab):2)
        {
          load.tab<-load.tab[order(abs(num.load.tab[,i-1]),decreasing=TRUE),]
          num.load.tab<-num.load.tab[order(abs(num.load.tab[,i-1]),decreasing=TRUE),]
        }
        if(gapped)  
        {
          load.tab<-
            do.call("rbind",
                    lapply(
                      1:nrow(load.tab),
                      function(i)
                      {
                        if(i %in% (1:nrow(num.load.tab))[apply(apply(num.load.tab==0,2,diff)==1,1,sum)==1])
                          rbind(load.tab[i,],"")
                        else load.tab[i,]
                      }))
        }
      }
      if(!is.null(bold.rows))
        load.tab[bold.rows,]<-gsub("<b></b>","",paste("<b>",load.tab[bold.rows,],"</b>",sep=""))
      factor.names<-gsub('Factor','Factor ',colnames(load.tab)[-1])
      tab.heads<-
        cbindTable(
          list(Text='Item',Other=list(align='left')),
          do.call("cbindTable",as.list(factor.names)))
      tab.rows<-
        do.call('cbindTable',
                lapply(
                  1:ncol(load.tab),
                  function(ii)list(Text=load.tab[,ii],Other=list(align=ifelse(ii==1,'left','center')))))
      notes<-
        paste(
          'Note: Loadings were estimated using the <em>',rotation,'</em> rotation. ',
          ifelse(threshold>0,paste('Only loadings exceeding &plusmn; ',threshold,' are shown.',sep=''),''),
          ifelse(
            is.numeric(boot),
            paste(' Confidence intervals were estimated from',boot-boot.quantiles[[3]],'bootstrap replicates.'),
            ''),
          sep='')
      tab.rows<-
        rbindTable(tab.rows,list(Text=notes,Other=list(align='left')))
      
      # make table of factor correlations
      
      if(!is.null(fa.out$rotmat)) 
      {
        tmat<-solve(fa.out$rotmat)
        fact.cor<-tmat %*% t(tmat)
        fact.cor<-cbind(factor.names,my.format(fact.cor,2))
        fact.cor.tab<-
          list(
            tab.heads=
              cbindTable(
                list(Text='',Other=list(align='left')),
                do.call("cbindTable",as.list(factor.names))),
            tab.rows=
              do.call('cbindTable',
                      lapply(
                        1:ncol(fact.cor),
                        function(j)list(Text=fact.cor[,j],Other=list(align=ifelse(j==1,'left','center'))))))                
      } else fact.cor.tab<-NULL
      list(tab.heads=tab.heads,tab.rows=tab.rows,fact.cor.tab=fact.cor.tab,load.tab=load.tab,fa.out=fa.out,n=n,na.pc=na.pc)
    }
  }





########################################################################################
#
# Functions to make figures
#


# function to duplicate scree.plot but with a neater simulation graphic - warning, it's slow if you increase simu!
#
# the difference is that scree.plot draws lines through eigenvalues from simu simulated random normal data sets 
# (same number of patients and variables as in the original data set, same pattern of missing data if use="P"),
# while neat.scree.plot permutes the original data set simu times, draws a solid line throught the means of the
# eigenvalues, and draws dotted lines at the +/- 95% confidence bounds. strictly speaking, this is not justified
# because the null distribution of each eigenvalue is dependent on all the others, but in practice it gives the same
# info as the scree.plot simulations but with a neater and more consistent output that doesn't obscure any data points.
# NB: simu=50 is fine for checking output but for a final report or paper it's worth increasing to >=500. 
#
# As of 2 May 2008 function also plots acceleration factor (2nd derivative) of the eigen values - peaks in this
# correspond to "accelerations" in the curve and the highest peak is the acceleration factor, which picks the "elbow"
# of the curve and is an objective solution to the Cattell subjective method.
#
# Changes to the function October 2011: Polychoric and polyserial correlations are used whenever >=1 of the 
# variables is an ordinal factor. Because these correlations are much slower to estimate than Pearson 
# correlations, the calculation of the null correlation band by permutation of the input data set became very 
# slow. I've therefore added an option that allows the null eigen values to be calculated from a data set of
# random independent normal variables. this is much faster but may underestimate the uncertainty 
# 

neat.scree.plot<-
  function(namefile,title="Scree Plot",simu=NULL,permute.cols=TRUE,accel=FALSE,method="pearson") 
  {
    if(any(sapply(namefile,is.ordered)))
    {
      require(polycor)
      print("Using polychoric correlations")
      cor.fun<-function(x)hetcor(x,std.err=FALSE,ML=FALSE,use="pairwise.complete",pd=TRUE)$correlations                        
    } else
    {
      cor.fun<-function(x)cor(x,use="pairwise.complete.obs",method=method)
    }
    ylab=ifelse(accel,"Eigenvalue\nAcceleration","Eigenvalue")
    eig.v<-eigen(cor.fun(namefile))$values
    acceleration<-if(accel) diff(diff(eig.v)) else NULL
    old.par<-par(mar=c(5.1,4.1,4.1,4.1))
    ylim<-range(c(acceleration,eig.v))
    plot(eig.v,type="n",xlab="Dimension",ylab="Eigenvalue",main="Scree plot",ylim=ylim)
    axis(1,at=1:length(eig.v),labels=FALSE,tcl=-0.25)
    if(is.numeric(simu)) 
    {   
      rand.ev<-
        do.call(
          'rbind',
          lapply(
            as.list(1:simu),
            function(i)
            {
              sim.cor<-
                if(permute.cols) cor.fun(permute(namefile)) else
                  cor(structure(rnorm(1:(ncol(namefile)*nrow(namefile))),dim=dim(namefile)))
              eigen(sim.cor,symmetric=TRUE,only.values=TRUE)$values
            }))
      q.rand.ev<-apply(rand.ev,2,function(x)quantile(x,c(0.025,0.975)))
      upper.q<-q.rand.ev[2,]
      lower.q<-q.rand.ev[1,]
      polygon(y=c(upper.q,rev(lower.q)),x=c(1:length(eig.v),length(eig.v):1),
              col=grey(0.7),border=NA)
      n.factors<-sum(eig.v>upper.q)
    } else n.factors<-NULL
    abline(h=1,lty=2)
    if(accel)
    {
      axis(4)  
      mtext("Acceleration",4,2,padj=1)
      points(x=1:length(acceleration)+1,y=acceleration,type="o",lty=3,pch=21,bg="white")
    }
    points(eig.v,type="o",pch=16)
    par(old.par)
    c(n.factors=n.factors)
  }


# function to permute columns of a data frame to make them independent. used by neat.scree.plot

permute<-
  function(data.list)
  {
    data.list<-as.data.frame(data.list)
    data.list<-
      lapply(data.list,sample,replace=FALSE)
    as.data.frame(data.list) 
  }




# function to remove html code from labels for plotting

plot.label<-
  function(x)
  {
    label(x)<-gsub("</sub>","",label(x))
    label(x)<-gsub("<sub>","",label(x))
    label(x)<-gsub("</sup>","",label(x))
    label(x)<-gsub("<sup>","",label(x))
    label(x)<-gsub("&ndash;","-",label(x))
    label(x)<-gsub("&gt;",">",label(x))
    label(x)<-gsub("&ge;",">=",label(x))
    label(x)<-gsub("&lt;","<",label(x))
    label(x)<-gsub("&le;","<=",label(x))          
    if(is.factor(x))
    {  
      levels(x)<-gsub("</sub>","",levels(x))
      levels(x)<-gsub("<sub>","",levels(x))
      levels(x)<-gsub("</sup>","",levels(x))
      levels(x)<-gsub("<sup>","",levels(x))
      levels(x)<-gsub("&ndash;","-",levels(x))
      levels(x)<-gsub("&gt;",">",levels(x))
      levels(x)<-gsub("&ge;",">=",levels(x))
      levels(x)<-gsub("&lt;","<",levels(x))
      levels(x)<-gsub("&le;","<=",levels(x))          
    }
    x
  }

# function to make file name and path for figure given figure number

fig.name<-
  function(fignum, figures.directory, fmt = ".jpg")
  {
    paste(
      figures.directory,'/figure',
      ifelse(regexpr('[.]',fignum)==2 | nchar(fignum)<2,paste('0',fignum,sep=''),fignum),
      fmt,sep='')
  }


# function to draw dotplot
# can adjust for covariates using the form: test=list(method="polr",adjust=dataset[,c("adj1","adj2")])
# CURRENTLY DOESN'T WORK FOR JUST ONE ADJUSTMENT COVARIATE E.G.test=list(method="polr",adjust=dataset[,"adj1"])
# e.g... 
#  dot.plot(x=av.data$adhd.inat,y=av.data$share.att.duration.tenths,x.factor=0.5,y.factor=1,
#    test=list(method="lm",adjust=av.data[,c("sex")]),pch=1,cex=0.6)
# ...won't adjust for sex

dot.plot<-
  function(x,y,x2=NULL,x.factor=1,y.factor=1,line.fn=median,
           test=list(method="wilcox",adjust=NULL),cex=1,cex.axis=1,pch=21,col="black",bg="white",...)
  {
    
    # function to turn 1,2,3,4,5,6,... into -1,1,-2,2,-3,3,...
    # used in splitting points to either side of central column 
    leftright<-function(x)ceiling(x/2)*((x/2==round(x/2))*2-1)
    
    # prepare data
    
    xy<-data.frame(x=x,y=y)
    xy$pch<-pch
    xy$col<-col
    xy$bg<-bg
    if(!is.null(test) && !is.null(test$adjust)) xy<-data.frame(xy,test$adjust)
    xy
    if(!is.null(x2)) xy$x2<-x2
    clean.xy<-na.omit(xy)
    clean.xy<-clean.xy[order(clean.xy$y),]
    x<-clean.xy$x
    y<-clean.xy$y
    x2<-clean.xy$x2
    if(!is.null(test) && !is.null(test$adjust)) test$adjust<-clean.xy[,names(test$adjust)]
    xplot<-
      if(is.null(x2)) x else
        factor(paste(x,x2,sep="\n"),apply(expand.grid(levels(x),levels(x2))[,1:2],1,paste,collapse="\n"))
    if(!is.null(x2)) label(xplot)<-""
    
    # prepare plot  
    old.par<-par(mar=c(6.1 - is.null(x2),4.1,4.1,2.1))
    ylim<-range(y)
    if(!is.null(test)) ylim[2]<-ylim[2]+diff(ylim)/10
    xlim<-range(as.numeric(xplot))+c(-0.5,0.5)
    plot(1,type="n",ylim=ylim,xlim=xlim,axes=FALSE,xlab=label(xplot),ylab=label(y),cex.lab=cex.axis,...)
    box()
    axis(1,at=sort(unique(as.numeric(xplot))),labels=levels(xplot),padj=0.5*(!is.null(x2)),
         las=ifelse(is.null(x2),0,2),cex.axis=cex.axis)
    axis(2,cex.axis=cex.axis)
    
    # define x and y gap sizes. I don't know if these will need much tuning for data of different scales
    ygap<-y.factor*diff(ylim)/50
    xgap<-x.factor*diff(xlim)/80
    
    # draw dots
    # the algorithm separates the sorted y points into groups of 0 to ygap, ygap to 2*ygap, 2*ygap to 3*ygap, etc.
    # the first point of each group is drawn centrally. subsequent points split to either side.
    
    
    by(data.frame(x=xplot,y=y,pch=clean.xy$pch,col=clean.xy$col,bg=clean.xy$bg),as.numeric(xplot),function(z)
    {
      z$x<-as.numeric(z$x)
      if(!is.null(line.fn)) lines(c(-0.3,0.3)+unique(z$x),rep(line.fn(z$y),2))
      i<-0
      counter<-0
      while(i<length(z$y))
      {
        i<-i+1
        j<-i
        counter<-counter+1
        while(j < length(z$y) & z$y[j+1]-z$y[i] < ygap)
        {
          j<-j+1
          z$x[j]<-z$x[j]+leftright(j-i)*xgap
        }
        i<-j 
      }
      points(x=z$x,y=z$y,cex=cex,pch=z$pch,col=z$col,bg=z$bg)
    })
    if(!is.null(test))
    {
      if(test$method=="wilcox") 
        p.value<-wilcox.test(y~x,exact=FALSE,correct=TRUE)$p.value
      if(test$method=="polr")
      {
        adj.text<-if(length(test$adjust)>0) paste("+test$adjust[[",1:length(test$adjust),"]]",sep="",collapse="") else NULL
        mod0<-polr(as.formula(paste("as.factor(y)~1",adj.text)))
        mod1<-polr(as.formula(paste("as.factor(y)~1+x",adj.text)))
        p.value<-anova(mod0,mod1)[2,"Pr(Chi)"]
      }
      if(test$method=="lm")
      {
        if(is.null(test$transform.fn)) test$transform.fn<-I
        adj.text<-if(length(test$adjust)>0) paste("+test$adjust[[",1:length(test$adjust),"]]",sep="",collapse="") else NULL
        #            adj.text<-if(length(test$adjust)>0) paste("+test$adjust",sep="",collapse="") else NULL
        
        print(test$adjust)
        print(adj.text)
        mod0<-lm(as.formula(paste("test$transform.fn(y)~1",adj.text)))
        mod1<-lm(as.formula(paste("test$transform.fn(y)~1+x",adj.text)))
        p.value<-anova(mod0,mod1)[2,"Pr(>F)"]
      }
      text(x=mean(xlim),y=ylim[c(2,2)]/1.05,paste(ifelse(is.null(x2),"",label(x)),p.format(p.value,p.eq=TRUE)))
    }
    par(old.par)
  }



# a function to fit a logistic regression model using elrm (exact logistic regression)
# this function does the pre-processing of the data.
# it also allows the MCMC chain to be extended until either the effective size
# (the number of MCMC samples adjusted to allow for autocorrelation) or the 
#  SE of the p-value estimate reaches some desired threshold.
# it won't work with interactions, unless they're coded as data columns.
# the output is the same as from elrm.
# most of the arguments are obvious or are explained at ?elrm. the more cryptic ones are:
# round.ndp: a list where each element name must be in names(source.data), and the elements
#   are the number of decimal places used to round the named variable. e.g. list(age = -1) 
#   will cause age to be rounded to the nearest decade. the reason for rounding is 
#   that continuous variables with mostly unique values tend to cause elrm to fail.
#   rounding is a crude way of reducing the number of unique values.
#   i find that rounding such that 6-8 unique values remain strikes a good balance between
#   achieving reliable results and preserving the information in the data.
#   the default is round.ndp=NULL.
# iter.used: the number of iterations after burnin in the first chain.
# iter.ext: MCMC can be extended to refine the estimates. iter.ext is the number 
#   of iterations in each extension chain.
# max.updates: max.updates is the maximum number of extensions allowed.
# max.p.value.se: MCMC will stop if the SE of the p.value falls below this threshold.
# min.effective.iter: MCMC will stop if the estimated effective number of iterations exceeds this threshold.
#   autocorrelation in the Markov chain will cause the effective number of iterations to be less
#   than iter.used. 
# example:
#   elrm.fit(outcome="adhd.inat",outcome.level="Yes",formula.rhs=~sex+factor1,interest=~factor1,
#     source.data=av.data,round.ndp=list(factor1=0),burnin=5000,iter.used=50000,iter.ext=10000,max.updates=30,
#     max.p.value.se=0.001,min.effective.iter=10000,r=2)
#


elrm.fit<-
  function(outcome,outcome.level,formula.rhs,interest,source.data,round.ndp=NULL,burnin=1000,
           iter.used=1000,iter.ext=1000,max.updates=0,min.effective.iter=0,max.p.value.se=Inf,r=4)
  {
    require(elrm)
    start.time<-Sys.time()
    
    # extract predictor terms from formula and make na-free data frame
    # with a binary 0/1 outcome
    
    x.terms<-attributes(terms(formula.rhs))$term.labels
    dat<-na.omit(source.data[,c(outcome,x.terms)])
    dat$y.label<-as.numeric(dat[,outcome]==outcome.level)
    dat[,outcome]<-NULL
    dat$n<-1
    
    # round predictors if necessary to avoid degeneracy
    
    for(n in names(round.ndp)) dat[,n]<-round(dat[,n],round.ndp[[n]])
    
    # collapse data frame to give binomial outcome
    
    binom.tab<-
      data.frame(
        do.call("rbind",
                by(dat,dat[,x.terms],
                   function(x)c(apply(x[,c("y.label","n")],2,sum),unlist(x[1,x.terms])))))
    names(binom.tab)[-1:-2]<-x.terms
    
    # fit model
    # needs to be constructed as an unevaluated expression, then evaluated
    # otherwise update.elrm won't run
    
    mod.formula<-formula(paste("y.label/n",paste(formula.rhs,collapse="")))
    mod.call<-
      substitute(
        elrm(mod.formula,interest=interest.formula,data=binom.tab,burnIn=burnin,iter=total.iter,r=r.par),
        list(mod.formula=mod.formula,burnin=burnin,total.iter=burnin+iter.used,interest.formula=interest,r.par=r))
    mod<-eval(mod.call)
    
    if(max.updates>0)
    {
      for(i in 1:max.updates)
      {
        print(c(EffectiveSize=effectiveSize(mod$mc),P=mod$p.values,P.SE=mod$p.values.se,autocorr(mod$mc)[,1,1]))        
        enough.samples <- mod$p.values.se<max.p.value.se | effectiveSize(mod$mc)>min.effective.iter
        if(enough.samples) break
        print(paste("Running Markov chain extension no",i))
        mod<-update(mod,iter=iter.ext)
      }
      print(paste("Number of Markov chain extensions:",i - 1 + !enough.samples))
      print(paste(ifelse(enough.samples,"At least one","None"),"of the criteria for p-value precision and effective MCMC size was met"))
    }  
    print(Sys.time()-start.time)
    mod
  }



# function to get kappa values from multiple raters. 
# I don't know how to deal with missing data for the moment so average over pairwise complete kappas
# if nan.rm=TRUE then uninformative pairwise comparisons which result in
# NaN (e.g. when there's no variation) are dropped to allow a mean kappa to be calculated.
# setting min.n causes pairwise kappa based on < min.n comparisons to be dropped.

kappa.multi<-
  function(rater.data,levels=NULL,nan.rm=FALSE,min.n=0)
  {
    require(vcd)
    if(is.null(levels(rater.data))) score.levels<-sort(unique(c(unlist(rater.data))))
    n.rater<-ncol(rater.data)
    mean.n.per.rater<-mean(apply(!is.na(rater.data),2,sum))
    k.list<-
      lapply(1:(n.rater-1),function(i)
      {
        k.tab<-
          sapply((i+1):n.rater,function(j)
          {
            rating.matrix<-
              table(
                factor(rater.data[,i],score.levels),
                factor(rater.data[,j],score.levels))
            n<-sum(rating.matrix)
            overall.kappa.stats<-Kappa(rating.matrix,weights="Fleiss-Cohen")
            overall.kappa<-overall.kappa.stats$Unweighted["value"]
            overall.kappa.quad<-overall.kappa.stats$Weighted["value"]
            if(n<min.n) overall.kappa.quad<-overall.kappa<-n<-NA
            c(k=overall.kappa,kw=overall.kappa.quad,mean.n=n)
          })
        t(k.tab)
      })
    c(apply(do.call("rbind",k.list),2,mean,na.rm=nan.rm),n.rater=n.rater,mean.n.per.rater=mean.n.per.rater)
  }




# funtion to adjust (approximately) logit-binomial GLMM predictions for jensen's inequality 
# method: mcculloch is more accurate, but the inverse function doesn't work, in which case use zeger
# p is the estimated probability, which should be biased where inverse = FALSE
# V is the sum of normal random effect variances
jensen.logit.adjust <-
  function(p, V, method = "mcculloch", inverse = FALSE) {
    stopifnot(!(method == "mcculloch" & inverse))
    Beta <- qlogis(p)
    if(method == "mcculloch") {
      return(plogis(Beta - 0.5 * V * tanh(Beta * (1 + 2 * exp(-0.5 * V))/6)))
    }
    if(method == "zeger") {
      if(inverse) 
        return(plogis(Beta * sqrt(256 * V / (75 * pi^2) + 1))) else
          return(plogis(Beta/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V)))
    }
  }


# Function to plot estimates with error bars
estimate.ci.plot <-
  function(x.term, y.term, plotdata, y.ci.lo = NULL, y.ci.hi = NULL, 
           x.at = NULL, x.labels = NULL, ylim = NULL, ...) {
    x <- plotdata[, x.term]
    y <- plotdata[, y.term]
    ci.lo <- plotdata[, y.ci.lo]
    ci.hi <- plotdata[, y.ci.hi]
    if(is.null(ylim) & !is.null(y.ci.lo) & !is.null(y.ci.hi)) 
      ylim <- c(min(ci.lo) * 0.95, max(ci.hi) * 1.05)
    plot(x, y, ylim = ylim, type = "n", axes = FALSE, ...)
    if(!is.null(y.ci.lo) & !is.null(y.ci.hi))
      arrows(x0 = x, x1 = x, y0 = ci.lo, y1 = ci.hi, angle = 90, code = 3, length = 0.1)
    points(x, y, ...)
    box()
    axis(1, at = x.at, labels = x.labels)
    axis(2)
  }


# Function to add square brackets around numbers
add_brackets <- function(labels) {
  gsub("(\\d)", "[\\1]", labels)
}


# Function to report messages (e.g. progress) from inside mclapply
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
