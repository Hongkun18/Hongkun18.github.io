library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
#library(shinyWidgets)
library(shinycssloaders)
library(shinythemes)
library(tidyverse)
library(shinyBS)

library(DT)
#library(readxl)
library(cowplot)
library(stringr)
library(plotly)

library(patchwork)
library(RColorBrewer)
library(randomcoloR)
#library(cowplot)
library(ggplot2)
library(devtools)

options(scipen=3,digits = 3)

#添加单细胞
library(Seurat)
library(ggsci)
library(viridis)

#
#drug_resistance dataset
#
pb<-readRDS("data/lung.rds")
table1 <- pb@meta.data[,c(1:6,36,38:41,46,54:56)]

#细胞通讯
library(ggplot2)

library(devtools)

# if (!require("BiocManager", quietly = TRUE))
   # install.packages("BiocManager")

# BiocManager::install("BiocGenerics")

# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")

# BiocManager::install("BiocNeighbors")

# # 要检测和安装的包名称列表
# packages_to_install <- c("car", "ragg", "sna", "ggpubr", "ggnetwork","NMF")

# # 检测并安装包
# for (package in packages_to_install) {
  # if (!require(package, quietly = TRUE)) {
    # install.packages(package)
  # }
# }

# packages_to_install2 <- c("presto", "circlize", "ComplexHeatmap", "CellChat")

# # 检测并安装包
# for (package in packages_to_install2) {
  # if (!require(package, quietly = TRUE)) {
    # install.packages(package)
  # }else {
  # # 如果已加载，直接使用
  # print("CellChat package is already loaded.")
	# }
# }
# 检查是否已加载 CellChat 包
# if (!require(CellChat, quietly = TRUE)) {
  # # 如果没有加载，尝试安装并加载
  # # install.packages('car')
  # # install.packages('ragg')
  # # install.packages('sna')
  # # install.packages('ggpubr')
  # # install.packages('ggnetwork')
  # # install.packages('NMF')
  # # devtools::install_github("immunogenomics/presto")
  # devtools::install_github("jokergoo/circlize")
  # devtools::install_github("jokergoo/ComplexHeatmap")
  # devtools::install_github("jinworks/CellChat")
  # library(CellChat)
	# } else {
  # # 如果已加载，直接使用
  # print("CellChat package is already loaded.")
# }

#library(CellChat)

# cellchat <- readRDS("data/cellchat_comparisonAnalysis_lungcancer_DR_vs_PR.rds")
# PR <- cellchat@LR$PR$LRsig
# PR$group <- "Partial_response"
# DR <- cellchat@LR$DR$LRsig
# DR$group <- "Drug_resistance"
# LR <- rbind(DR,PR)
# rownames(LR) <- NULL
# #将pb与LR放在一起
# pb@assays$RNA@meta.features<- LR

#
#NSCLC dataset
#
#TNM_I
#nsclc<- readRDS("data/1.paired_sample_core_atlas_luad_TNM_I_corrected_clinical_information_smoking_seurat2.rds")
#TNM-II
# nsclc<- readRDS("data/1.paired_sample_core_atlas_luad_TNM_II_corrected_clinical_information_smoking_seurat2.rds")
# nsclc$celltype.main <- nsclc$ann_coarse
# nsclc <- SetIdent(nsclc, value = "celltype.main")
# nsclc$sample <- as.character(nsclc$sample)

## table
load("data/Table.RData")
#table1[,-c(1:2,5)] <- signif(table1[,-c(1:2,5)],3)
table2A[,5:7] <- signif(table2A[,5:7],3)

table2B[,5:9] <- signif(table2B[,5:9],3)
table2C[,6:10] <- signif(table2C[,6:10],3)
table3A[,9:10] <- signif(table3A[,9:10],3)
table3B[,10:14] <- signif(table3B[,10:14],3)

#single cell 

#
#drug resistance dataset
#
pb<-readRDS("data/lung.rds")
IdentRename<-function(pb,oldname,newname){
  cluster.ids <-levels(pb)
  cluster.ids[which(cluster.ids  == oldname)] <- newname
  names(cluster.ids) <- levels(pb)
  pb <- RenameIdents(pb, cluster.ids)
  #pbmc<-IdentRename(pb,oldname,newname)
}


# Home 模块
tab_home <- dashboardPage(
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  dashboardBody(
     fluidRow(
       box(#title = tags$h3("HMOX1 Expression Influence the Role of Macrophage in EGFR-TKI Resistance of Lung Adenocarcinoma"),
           solidHeader = F, collapsible = F, width = 12, 
           tabPanel(title=NULL,
					box(title = tags$p("LUAD", style = "margin-top: 5px; margin-left: 10px;font-size: 26px;"),
                        solidHeader = TRUE, collapsible = FALSE, width = 6,
                        #tags$p(img(src = "img/workflow1.png", style='width:90%;height:90%;align:center')),
						tags$p(img(src = "img/draft_figure_final_version.png", style='width:100%;height:100%;text-align:center')),
                          style='border:2px solid #C8C8C8;border-radius: 10px;text-align:justify;'),
                    #box(title = tags$p("Workflow", style = "margin-top: 5px; margin-left: 10px;font-size: 26px;"),
					box(title = tags$p("NSCLC", style = "margin-top: 5px; margin-left: 10px;font-size: 26px;"),
                        solidHeader = TRUE, collapsible = FALSE, width = 6,
                        #tags$p(img(src = "img/workflow1.png", style='width:90%;height:90%;align:center')),
						tags$p(img(src = "img/NSCLC.png", style='width:100%;height:100%;text-align:center')),
                          style='border:2px solid #C8C8C8;border-radius: 10px;text-align:justify;')
                  ),
           h3(strong(".")),
		   tags$p('The dataset from the following publication:', style = "font-size: 150%;"),
           #tags$p('Please cite the following publication:', style = "font-size: 150%;"),
           tags$p(a('1.Salcher S, et al. High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer[J]. Cancer cell, 2022, 40(12): 1503-1520.e8.',href="https://www.cell.com/cancer-cell/pdf/S1535-6108(22)00499-8.pdf", target="view_window")),
		   tags$p(a('2.Hongkun Fang et al. HMOX1 Expression Influence the Role of Macrophage in EGFR-TKI Resistance of Lung Adenocarcinoma.',#Cancer Research.2021 Aug 15;81(16):4205-4217.
		   href="https://cancerres.aacrjournals.org/content/81/16/4205", target="view_window")
                  )
           )
       )
  )
)

#Document 模块 
document <- dashboardPage(
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  dashboardBody(
    box(
    title = NULL, solidHeader = TRUE, collapsible = FALSE,
    width = 12,
    tags$ul(
      htmlOutput("result")
    )
    )
  )
)

#Contact 模块
contact <- dashboardPage(
  
  dashboardHeader(disable = T), 
  dashboardSidebar(disable = T), 
  dashboardBody(
    fluidRow(
    box(
      title = NULL, solidHeader = TRUE, collapsible = FALSE,
      width = 12, # solidHeader=TRUE can remove the top boarder
      tags$h1(strong("Contact us")),
      
      tags$p(strong("Please don't hesitate to address comments, questions or suggestions regarding this website.", 
             style = "font-size: 150%;")),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),
      tags$h2("Corresponding author"),
      tags$h4("Qiyuan Li: qiyuan.li@xmu.edu.cn"),
	  tags$h4("Ying Zhou: yingzhou@xmu.edu.cn"),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),
      tags$h2("First author"),
      tags$h4("Hongkun Fang: hongkunfang@126.com"),
      tags$hr(style="border-top: 1px dashed #A9A9A9"),
      tags$h2("Institute"),
      tags$h4(a("National Institute for Data Science in Health and Medicine, Xiamen University", href="https://nidshm.xmu.edu.cn/",target="view_window")),
      tags$h4(a("School of Medicine Xiamen University, China", href="https://med.xmu.edu.cn/", target="view_window")),
      tags$hr(style="border-top: 1px dashed #A9A9A9")
    ),
    box(
      title = NULL, solidHeader = TRUE, collapsible = FALSE,
      width = 12, # solidHeader=TRUE can remove the top boarder
      dashboardFooter(
        left = HTML('<footer><h10>visits since Nov 30, 2020</h10>
                          <script type="text/javascript" src="//rf.revolvermaps.com/0/0/6.js?i=5btpst0tngq&amp;m=7&amp;c=e63100&amp;cr1=ffffff&amp;f=arial&amp;l=0&amp;bv=90&amp;lx=-420&amp;ly=420&amp;hi=20&amp;he=7&amp;hc=a8ddff&amp;rs=80" async="async"></script>
                          </footer>')
        )
      )
    )
  )
)

###### UI start
ui <- fluidPage(
  useShinyjs(),
  #extendShinyjs(text = jscode, functions = "refresh"),
  tags$head(tags$meta(name = "viewport", content = "width=1260")),
  tags$style(HTML("
          .navbar .navbar-nav {float: right}
          .navbar .navbar-header {float: right}
          .navbar {background-color:#2C3E50;}
          .navbar {font-color:black;}		

		  # .my-image {
		  # max-width: 100%;  /* 设置最大宽度为父容器的100% */
		  # max-height: 100%; 
		  # height: auto;       /* 根据高度等比例调整宽度 */
		  # width: 100%;       /* 设置宽度为高度的50% */
          # display: block;   /* 将图像设置为块级元素，以便可以设置上下边距 */
          # margin: 0 auto;   /* 居中显示 */}
        "
		)),
  navbarPage(
    #title = NULL,
	title = "Single Cell Lung Cancer",
    id = 'navbar',
    #windowTitle = "Genetic Determinants of the Somatic Selection of Mutational Processes in 3,566 Human Cancers",
    theme = shinytheme("flatly"),
    tabPanel("Home", value = 'home', tab_home,icon=icon('home')),
    
	
	#Explorer 模块，添加多个模块
	tabPanel("Explorer", 
	         icon = icon('chart-area'),  # 设置图标样式为 "chart-area"（这是图形按钮的一种示例）
	         #actionButton("exploreButton", "", icon = icon("search"), width = "100px", height = "100px"),
	        # verbatimTextOutput("explorationResult"),
			 
	#添加单细胞展示网页
	#shinyUI(fluidPage(
  
	  # Application title
	  titlePanel(p("SeuratReport" , style = "color:#3474A7")),
	  
	  # Sidebar with a slider input for number of bins 
	  sidebarLayout(
		 		
		sidebarPanel(
		  
			##################
		  fileInput('file1', 'Choose Seurat  object  RDS File',
					accept=c('rds', 
							 '', 
							 '.rds')),
		  conditionalPanel(
					condition = "input.sample == ture",
					selectInput("ChooseSample",'ChooseSample',as.list( c("All",unique(c(pb$sample))))) #,nsclc$sample))))) #  levels(pb)  as.list(pb@meta.data$seurat_clusters)

		 ),
		  #添加umap tsne显示方式
		  # conditionalPanel(
			# condition = "input.smoother == ture",
			# selectInput("comment","reduction:",list("umap_harmony","umap","tsne","pca"))
		  # ),
		  
		  # textInput(inputId = "clustername",
					# label =  "clusterName",
					# value = "HERE"
		  # ),

		  # conditionalPanel(
					# condition = "input.cluster == ture",
					# selectInput("Choosecluster",'ChooseCluster',as.list( c("All",levels(pb) ))) #  levels(pb)  as.list(pb@meta.data$seurat_clusters)

		  # ),
		  
		  conditionalPanel(
			condition = "input == ture",
			selectInput("gene",'ChooseGene',as.list(c(rownames(pb))))#,rownames(nsclc@assays$RNA@counts))))# unique(c(rownames(pb),rownames(nsclc)))
			
		  ),
		 
		  conditionalPanel(
			condition = "input == ture",
			#selectInput("gene",'ChooseGene',as.list( rownames(pb) ))
			textAreaInput("geneinput", "Dotplot Gene", "")
			
		  ),
		  
		  
		  
		  #添加保存图片的方法
		  #br(), #插入一个换行
		  hr(), #插入一条水平线
		  conditionalPanel(
					condition = "input.sample == ture",
					selectInput("plotType",'Select Plot Type',c("DimPlot","VlnPlot","FeaturePlot", "DotPlot","Cellchat","Cellchat2"))),
		  numericInput("plotWidth", "Plot Width", 800),
		  numericInput("plotHeight", "Plot Height", 400),
		  selectInput("fileType", "Save As", c("png", "jpeg", "pdf"), selected = "png"),
		  downloadButton("downloadPlot", "Download Plot"),
		  #br(), #插入一个换行
		  hr(), #插入一条水平线
		  actionButton("updatePlot", "Update Plot"),
		  
		  ),
		
		#主面板
		# Show a plot of the generated distribution
		mainPanel(
		  #h3("DimPlot"),
		  tabsetPanel(type = "tabs",
			#LUAD
			tabPanel("LUAD",
				tableOutput('contents'),
				imageOutput("p1", width = "auto", height = "auto", click = NULL, #width = "100%",height = "400px"
							dblclick = NULL, hover = NULL,
							brush = NULL, inline = FALSE),
				tabsetPanel(
					tabPanel("Overview",plotOutput("DimPlot")),
					tabPanel("VlnPlot",plotOutput("VlnPlot")),
					tabPanel("FeaturePlot",plotOutput("FeaturePlot")),
					tabPanel("DotPlot",plotOutput("DotPlot")),
					#图片代替绘图
					#tabPanel("Cell-Cell Communication",plotOutput("Cellchat")),					
					tabPanel("Cell-Cell Communication",
						tags$p(img(src = "img/LUAD_cellchat1.png", style='width:100%;height:100%;text-align:center')),
                          style='border:0px solid #C8C8C8;border-radius: 10px;text-align:justify;'),
														
					#tabPanel("Cell-Cell LR Pathway",plotOutput("Cellchat2",height = "1200px",width = "1000px",inline = FALSE))
					tabPanel("Cell-Cell LR Pathway",
						tags$p(img(src = "img/LUAD_cellchat2.png", style='width:100%;height:100%;text-align:center')),
                          style='border:0px solid #C8C8C8;border-radius: 10px;text-align:justify;'),
					
				)
		  ),
			# #NSCLC
			# tabPanel("NSCLC",
				# tableOutput('contents1'),
				# imageOutput("p1_nsclc", width = "auto", height = "auto", click = NULL, #width = "100%",height = "400px"
							# dblclick = NULL, hover = NULL,
							# brush = NULL, inline = FALSE),
				# tabsetPanel(
					# tabPanel("Overview",plotOutput("DimPlot1")),
					# tabPanel("VlnPlot",plotOutput("VlnPlot1")),
					# tabPanel("FeaturePlot",plotOutput("FeaturePlot1")),
					# tabPanel("DotPlot",plotOutput("DotPlot1")),
					# #图片代替绘图
					# #tabPanel("Cell-Cell Communication",plotOutput("Cellchat3")),
					# #tabPanel("Cell-Cell LR Pathway",plotOutput("Cellchat4",height = "1200px",width = "1000px",inline = FALSE))
									
					# tabPanel("Cell-Cell Communication",
						# tags$p(img(src = "img/NSCLC_cellchat1.png", style='width:100%;height:100%;text-align:center')),
                          # style='border:0px solid #C8C8C8;border-radius: 10px;text-align:justify;'),
														
					# tabPanel("Cell-Cell LR Pathway",
						# tags$p(img(src = "img/NSCLC_cellchat2.png", style='width:100%;height:100%;text-align:center')),
                          # style='border:0px solid #C8C8C8;border-radius: 10px;text-align:justify;'),
					
					
					
					
				# )
			# )
		)
	  )			 		 
	 )

	),
	
	
	#MetaData 添加LUAD和NSCLC模块
	 tabPanel("MetaData",value ="table1",icon = icon('table'),
		 tabsetPanel(type = "tabs",
                  tabPanel("LUAD", value = "table1",
					tabsetPanel(type = "tabs",
						tabPanel("metadata",value = "table1",
							column(4,
							selectInput("sample1",
										"Sample:",
										c("All",
										  unique(as.character(pb$sample))))
							),
							column(4,
								selectInput("celltype1",
											"Celltype:",
											c("All",
											  unique(as.character(pb$celltype.main))))
							),
							column(4,
								selectInput("Group1",
											"Group:",
											c("All",
											  unique(as.character(pb$Group))))
							),
							   DT::dataTableOutput("table1"),
							   style = "overflow-x: scroll;"
					  ),
						tabPanel("Significant niLigand-Receptor Interactioons Pairs",value = "table2",
							column(4,
							selectInput("pathway1",
										"Pathway:",
										c("All",
										  unique(as.character(pb@assays$RNA@meta.features$pathway_name))))
							),
							column(4,
								selectInput("Group2",
											"Group:",
											c("All",
											  unique(as.character(pb@assays$RNA@meta.features$group))))
							),
                           DT::dataTableOutput("table2"),
                           style = "overflow-x: scroll;"
						 )
					  )
					),
						
				# tabPanel("NSCLC", value = "table3",
									# tabsetPanel(type = "tabs",
										# tabPanel("metadata",value = "table3",
											# column(4,
											# selectInput("sample1",
														# "Sample:",
														# c("All",
														  # unique(as.character(nsclc$sample))))
											# ),
											# column(4,
												# selectInput("celltype1",
															# "Celltype:",
															# c("All",
															  # unique(as.character(nsclc$celltype.main))))
											# ),
											# column(4,
												# selectInput("Group1",
															# "Group:",
															# c("All",
															  # unique(as.character(nsclc$origin))))
											# ),
											   # DT::dataTableOutput("table3"),
											   # style = "overflow-x: scroll;"
									  # ),
										# tabPanel("Significant niLigand-Receptor Interactioons Pairs",value = "table4",
											# column(4,
											# selectInput("pathway1",
														# "Pathway:",
														# c("All",
														  # unique(as.character(nsclc@assays$RNA@meta.features$pathway_name))))
											# ),
											# column(4,
												# selectInput("Group2",
															# "Group:",
															# c("All",
															  # unique(as.character(nsclc@assays$RNA@meta.features$group))))
											# ),
										   # DT::dataTableOutput("table4"),
										   # style = "overflow-x: scroll;"
										 # )
									  # )
									# ),
			)
	 
	 ),


	
	#Download 模块
    tabPanel("Download", icon = icon('download'),## icon:设置图标样式为 "download"
		titlePanel("Downloading Data"),
		column(12, style = "text-align: left;"),
		# tags$style(type="text/css", 
                      # ".tab-content {height: 70vh;width: 70%;}"),  # 设置tabPanel的高度和宽度			  
		sidebarLayout(

			# Sidebar panel for inputs ----
			sidebarPanel(

			  # Input: Choose dataset ----
			  selectInput("dataset", "Choose a dataset:",
						  #choices = c("rock", "pressure", "cars")),
						  choices = c("LUAD")#, "NSCLC")
						  ),
			  # Input: Choose dataset  datatype----			  
			  selectInput("datatype", "Choose a datatype:",
						  #choices = c("rock", "pressure", "cars")),
						  # choices = c("Meta data", "Expression Matrix","Significant Ligand-Receptor Interaction Pairs")
						  choices = c("Meta data", "Significant Ligand-Receptor Interaction Pairs")
						  ),
				# Button  
				#输出形式选择
			  radioButtons("filetype", "File type:",
                   choices = c("txt","csv")),	
			  downloadButton("downloadData", "Download")  ##downloadButton("Download", "Download.csv")), 

			),
			# Main panel for displaying outputs ----
			mainPanel(
			  #tableOutput("table") #table的展示方式
			  DT::dataTableOutput("table"),
			  style = "overflow-x: auto;"
				)
			),	
	),
	#Contact 模块
	tabPanel("Contact",value = 'contact', contact,icon = icon("envelope"))
	),
	# Add custom CSS to move the title to the left
	tags$head(
	  tags$style(HTML("
		.navbar-header {float: left !important; margin-left: 15px !important;}
		.navbar-brand {padding: 15px 0 0 60px !important;} #title往右移动60px
	  "))
	)
	
)

# shinyUI(ui)
###### UI end

server <- function(input, output, session) { 
  seed <- reactiveVal()
  
  ######单细胞数据
  #
  #1.drug resistance
  #
  output$contents <-  DT::renderDataTable({
		inFile <- input$file1
		if (is.null(inFile)){
		  return(NULL)
		}else{
		  pb<- readRDS(inFile$datapath)
		  pb
		  req(pb)
		  
		}
	  })
	   
	  #Overview DimPlot
	  colourCount = length(unique(pb$celltype.main))
	  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	  
	  output$DimPlot <- renderPlot({
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster) #input$clustername
		if(input$ChooseSample=="All"){
		#p1 <- DimPlot(c, label = FALSE,reduction   = input$comment) #+ NoLegend()

		p1 = DimPlot(c, group.by="celltype.main",reduction="umap_harmony",pt.size = 2)+ #reduction=input$comment
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				scale_color_manual(values = getPalette(colourCount))#+NoLegend()
	
		#p2 <- DimPlot(c, label = FALSE,reduction   = input$comment,group.by="Group")
		p2 <- DimPlot(c,label = FALSE,group.by = "Group",reduction="umap_harmony", cols = c("#1b9e77","#d95f02"), pt.size =1)+
        labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		print(p1+p2)
		#DimPlot(c, label = TRUE,reduction   = input$comment) + NoLegend()
		}else{
		d <- subset(c,sample==input$ChooseSample)
		p1 <- DimPlot(d, group.by="celltype.main",reduction="umap_harmony",pt.size = 2)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				scale_color_manual(values = getPalette(colourCount))
		#print(p1)
		#选择颜色
		if(unique(d$Group)=="Drug_resistance"){color="#1b9e77"}else{color="#d95f02"}
		p2 <- DimPlot(d,label = FALSE,group.by = "Group",reduction="umap_harmony", cols = color, pt.size =1)+
        labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		#p2 <- DimPlot(d, label = FALSE,reduction   = input$comment,group.by="Group")
		print(p1+p2)
		}
		#DimPlot(IdentRename(pb,input$Choosecluster,input$clustername), label = TRUE,reduction   = input$comment) + NoLegend()

	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  
	  #VlnPlot
	  output$VlnPlot<-renderPlot({
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
		if(input$ChooseSample=="All"){
		VlnPlot(c, features = input$gene, pt.size = 0.2, ncol = 1)+scale_fill_manual(values = getPalette(colourCount))
		}else{
		d <- subset(c,sample==input$ChooseSample)
		VlnPlot(d, features = input$gene, pt.size = 0.2, ncol = 1)+scale_fill_manual(values = getPalette(colourCount))
		}
		#VlnPlot(pb, features = input$gene, pt.size = 0.2, ncol = 1)
	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  

	  #FeaturePlot	  
	  output$FeaturePlot<-renderPlot({
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
		if(input$ChooseSample=="All"){
		FeaturePlot(c, features = input$gene, reduction="umap_harmony", pt.size = 0.2, ncol = 1)
		}else{
		d <- subset(c,sample==input$ChooseSample)
		FeaturePlot(d, features = input$gene, reduction="umap_harmony", pt.size = 0.2, ncol = 1)
		}
		#FeaturePlot(IdentRename(pb,input$Choosecluster,input$clustername), features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  
	  #DotPlot
	  output$DotPlot<-renderPlot({
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
		gene_names <- strsplit(input$geneinput, "\n")[[1]]
		if(input$ChooseSample=="All"){
		DotPlot(c, features = rev(as.character(unique(gene_names))))+ 
		theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        ylab("Cell type")+
        xlab("Marker gene")+
        coord_flip()+
        scale_color_viridis()
		}else{
		d <- subset(c,sample==input$ChooseSample)
		DotPlot(d, features = rev(as.character(unique(gene_names))))+ 
		theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        ylab("Cell type")+
        xlab("Marker gene")+
        coord_flip()+
        scale_color_viridis()
		}
	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  
	  #------CellChat文件
	  #Overview DimPlot
	  #colourCount = length(unique(pb$celltype.main))
	  #getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	  #Cell-Cell communication
	  output$Cellchat<-renderPlot({
		
		#2.1.1数量与强度差异网络图
		par(mfrow = c(1,2), xpd=TRUE)
		netVisual_diffInteraction(cellchat, weight.scale = T)
		netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
		
		#2.1.2数量与强度差异热图
		#gg1 <- netVisual_heatmap(cellchat)
		#gg2 <- netVisual_heatmap(cellchat, measure = "weight")
		#> Do heatmap based on a merged object
		#gg1 + gg2	
	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  #Cell-Cell Pathway
	  output$Cellchat2<-renderPlot({
		# c <- IdentRename(cellchat,input$Choosecluster,input$Choosecluster)
		# if(input$ChooseSample=="All"){
		# FeaturePlot(c, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# }else{
		# d <- subset(c,sample==input$ChooseSample)
		# FeaturePlot(d, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# }
		
		#LRsig <- setdiff(cellchat@LR$PR$LRsig$pathway_name,cellchat@LR$DR$LRsig$pathway_name)
		LRsig <-intersect(cellchat@LR$PR$LRsig$pathway_name,cellchat@LR$DR$LRsig$pathway_name)
		LRsig_top10 <- LRsig[1:10]
		netVisual_bubble(cellchat, sources.use = c(2,3,6), targets.use = c(1,4,5,7,8,9),signaling = LRsig_top10,  comparison = c(1, 2), angle.x = 45)
		# gg1 <- netVisual_bubble(cellchat, sources.use = c(2,3,6), targets.use = c(1,4,5,7,8,9),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PR", angle.x = 45, remove.isolate = T)
		# gg2 <- netVisual_bubble(cellchat, sources.use =c(2,3,6), targets.use =  c(1,4,5,7,8,9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PR", angle.x = 45, remove.isolate = T)
		# gg1+gg2
		
	  },
	  #输出图片的长宽
	  width = function() input$plotWidth, height = function() input$plotHeight
	  )
	  
	  
	  
  #
  #2.NSCLC
  #
  # output$contents1 <-  DT::renderDataTable({
		# inFile <- input$file1
		# if (is.null(inFile)){
		  # return(NULL)
		# }else{
		  # nsclc <- readRDS(inFile$datapath)
		  # nsclc
		  # req(nsclc) 
		# }
	  # })
	  
	  # #Overview DimPlot
	  # library(randomcoloR)
	  # colourCount1 = length(unique(nsclc$celltype.main))
	  # palette <- distinctColorPalette(colourCount1)
	  
	  # output$DimPlot1 <- renderPlot({
		# c <- IdentRename(nsclc,input$Choosecluster,input$Choosecluster) #input$clustername
		# if(input$ChooseSample=="All"){
		# #p1 <- DimPlot(c, label = FALSE,reduction   = input$comment) #+ NoLegend()

		# p1 = DimPlot(c, group.by="celltype.main",pt.size = 2, label=F, raster=FALSE)+ #reduction=input$comment,
				# labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				# scale_color_manual(values = palette)
				# #scale_color_manual(values = getPalette(colourCount))#+NoLegend()
	
		# #p2 <- DimPlot(c, label = FALSE,reduction   = input$comment,group.by="Group")
		# p2 <- DimPlot(c,label = FALSE,group.by = "origin", cols = c("#1b9e77","#d95f02"), pt.size =1, raster=FALSE)+ #reduction=input$comment,
        # labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		# print(p1+p2)
		# #DimPlot(c, label = TRUE,reduction   = input$comment) + NoLegend()
		# }else{
		# d <- subset(c,sample==input$ChooseSample)
		# p1 <- DimPlot(d, group.by="celltype.main",pt.size = 2, raster=FALSE)+#reduction=input$comment,
				# labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				# scale_color_manual(values = palette)
				# #scale_color_manual(values = getPalette(colourCount))
		# #print(p1)
		# #选择颜色
		# if(unique(d$origin)=="normal_adjacent"){color="#1b9e77"}else{color="#d95f02"}
		# p2 <- DimPlot(d,label = FALSE,group.by = "origin",reduction=input$comment, cols = color, pt.size =1)+
        # labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		# #p2 <- DimPlot(d, label = FALSE,reduction   = input$comment,group.by="Group")
		# print(p1+p2)
		# }
		# #DimPlot(IdentRename(pb,input$Choosecluster,input$clustername), label = TRUE,reduction   = input$comment) + NoLegend()

	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  
	  # #VlnPlot
	  # output$VlnPlot1<-renderPlot({
		# c <- IdentRename(nsclc,input$Choosecluster,input$Choosecluster)
		# if(input$ChooseSample=="All"){
		# VlnPlot(c, features = input$gene, pt.size = 0.2, ncol = 1)+scale_fill_manual(values = palette)
		# }else{
		# d <- subset(c,sample==input$ChooseSample)
		# VlnPlot(d, features = input$gene, pt.size = 0.2, ncol = 1)+scale_fill_manual(values = palette)
		# }
		# #VlnPlot(pb, features = input$gene, pt.size = 0.2, ncol = 1)
	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  

	  # #FeaturePlot	  
	  # output$FeaturePlot1<-renderPlot({
		# c <- IdentRename(nsclc,input$Choosecluster,input$Choosecluster)
		# if(input$ChooseSample=="All"){
		# FeaturePlot(c, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# }else{
		# d <- subset(c,sample==input$ChooseSample)
		# FeaturePlot(d, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# }
		# #FeaturePlot(IdentRename(pb,input$Choosecluster,input$clustername), features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  
	  # #DotPlot
	  # output$DotPlot1<-renderPlot({
		# c <- IdentRename(nsclc,input$Choosecluster,input$Choosecluster)
		# gene_names <- strsplit(input$geneinput, "\n")[[1]]
		# if(input$ChooseSample=="All"){
		# DotPlot(c, features = rev(as.character(unique(gene_names))))+ 
		# theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        # ylab("Cell type")+
        # xlab("Marker gene")+
        # coord_flip()+
        # scale_color_viridis()
		# }else{
		# d <- subset(c,sample==input$ChooseSample)
		# DotPlot(d, features = rev(as.character(unique(gene_names))))+ 
		# theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        # ylab("Cell type")+
        # xlab("Marker gene")+
        # coord_flip()+
        # scale_color_viridis()
		# }
	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  
	  
	  # #------CellChat文件
	  # #Overview DimPlot
	  # #Cell-Cell communication
	  # output$Cellchat3<-renderPlot({
		
		# #2.1.1数量与强度差异网络图
		# par(mfrow = c(1,2), xpd=TRUE)
		# netVisual_diffInteraction(cellchat_nsclc, weight.scale = T)
		# netVisual_diffInteraction(cellchat_nsclc, weight.scale = T, measure = "weight")
		
		# #2.1.2数量与强度差异热图
		# #gg1 <- netVisual_heatmap(cellchat)
		# #gg2 <- netVisual_heatmap(cellchat, measure = "weight")
		# #> Do heatmap based on a merged object
		# #gg1 + gg2
	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  
	  # #Cell-Cell Pathway
	  # output$Cellchat4<-renderPlot({
		# # c <- IdentRename(cellchat,input$Choosecluster,input$Choosecluster)
		# # if(input$ChooseSample=="All"){
		# # FeaturePlot(c, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# # }else{
		# # d <- subset(c,sample==input$ChooseSample)
		# # FeaturePlot(d, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		# # }
		
		# #LRsig <- setdiff(cellchat@LR$PR$LRsig$pathway_name,cellchat@LR$DR$LRsig$pathway_name)
		# LRsig <-intersect(cellchat_nsclc@LR$NAT$LRsig$pathway_name,cellchat_nsclc@LR$TNM_II$LRsig$pathway_name)
		# LRsig_top100 <- LRsig[1:100]
		# netVisual_bubble(cellchat_nsclc, sources.use = c(1,4,6,9), targets.use = c(2,3,5,7,8),signaling = LRsig_top100,  comparison = c(1, 2), angle.x = 45)
		# # gg1 <- netVisual_bubble(cellchat, sources.use = c(2,3,6), targets.use = c(1,4,5,7,8,9),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in PR", angle.x = 45, remove.isolate = T)
		# # gg2 <- netVisual_bubble(cellchat, sources.use =c(2,3,6), targets.use =  c(1,4,5,7,8,9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in PR", angle.x = 45, remove.isolate = T)
		# # gg1+gg2
		
	  # },
	  # #输出图片的长宽
	  # width = function() input$plotWidth, height = function() input$plotHeight
	  # )
	  
	
	#
	#
	#
	#
	
	  
	# 下载plot按钮的功能
	  output$downloadPlot <- downloadHandler(
		filename = function() {
		  paste("",input$plotType,".", input$fileType, sep = "")
		},
		content = function(file) {
		#添加类型
		if(input$fileType == 'pdf'){
			pdf(file,width =input$plotWidth/100, height =input$plotHeight/100)
		  }else if(input$fileType == 'png'){
			png(file,width =input$plotWidth, height =input$plotHeight)
		  }else{
			jpeg(file,width =input$plotWidth, height =input$plotHeight)
		  }
		#plot
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster) #input$clustername
		if(input$ChooseSample=="All"){
			if (input$plotType == "FeaturePlot") {
				p <- FeaturePlot(c, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
			}else if(input$plotType == "DotPlot"){
				gene_names <- strsplit(input$geneinput, "\n")[[1]]
				p <- DotPlot(c, features = rev(as.character(unique(gene_names))))+ 
				theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
				ylab("Cell type")+
				xlab("Marker gene")+
				coord_flip()+
				scale_color_viridis()
				}else if(input$plotType == "VlnPlot"){
				p <- VlnPlot(c, features = input$gene, pt.size = 0.2, ncol = 1)
				}else if(input$plotType == "DimPlot"){
				#p <-DimPlot(c, label = FALSE,reduction   = input$comment)
				p1 = DimPlot(c, group.by="celltype.main",reduction=input$comment,pt.size = 2)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				scale_color_manual(values = getPalette(colourCount))#+NoLegend()

				p2 <- DimPlot(c,label = FALSE,group.by = "Group",reduction=input$comment, cols = c("#1b9e77","#d95f02"), pt.size =1)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
				p <- p1+p2
				
				}

		}else{
		d <- subset(c,sample==input$ChooseSample)
		if (input$plotType == "FeaturePlot") {
				p <- FeaturePlot(d, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
			}else if(input$plotType == "DotPlot"){
				gene_names <- strsplit(input$geneinput, "\n")[[1]]
				p <- DotPlot(d, features = rev(as.character(unique(gene_names))))+ 
				theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
				ylab("Cell type")+
				xlab("Marker gene")+
				coord_flip()+
				scale_color_viridis()
				}else if(input$plotType == "VlnPlot"){
				p <- VlnPlot(d, features = input$gene, pt.size = 0.2, ncol = 1)
				}else if(input$plotType == "DimPlot"){
				#p <-DimPlot(d, label = FALSE,reduction   = input$comment)
				p1 <- DimPlot(d, group.by="celltype.main",reduction=input$comment,pt.size = 2)+
						labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
						scale_color_manual(values = getPalette(colourCount))
				#print(p1)
				#选择颜色
				if(unique(d$Group)=="Drug_resistance"){color="#1b9e77"}else{color="#d95f02"}
				p2 <- DimPlot(d,label = FALSE,group.by = "Group",reduction=input$comment, cols = color, pt.size =1)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
				#p2 <- DimPlot(d, label = FALSE,reduction   = input$comment,group.by="Group")
				p<- p1+p2
				
				}
		
		}
		print(p)
		dev.off()
		  
		}
	  )
	  
	  
	  # 更新plot的按钮功能
	  observeEvent(input$updatePlot, {
	    #overview Dimplot
		output$DimPlot <- renderPlot({
		c <- IdentRename(pb,input$Choosecluster,input$Choosecluster) #input$clustername
		if(input$ChooseSample=="All"){
		#p1 <- DimPlot(c, label = FALSE,reduction   = input$comment) #+ NoLegend()

		p1 = DimPlot(c, group.by="celltype.main",reduction=input$comment,pt.size = 2)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				scale_color_manual(values = getPalette(colourCount))#+NoLegend()
	
		#p2 <- DimPlot(c, label = FALSE,reduction   = input$comment,group.by="Group")
		p2 <- DimPlot(c,label = FALSE,group.by = "Group",reduction=input$comment, cols = c("#1b9e77","#d95f02"), pt.size =1)+
        labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		print(p1+p2)
		#DimPlot(c, label = TRUE,reduction   = input$comment) + NoLegend()
		}else{
		d <- subset(c,sample==input$ChooseSample)
		p1 <- DimPlot(d, group.by="celltype.main",reduction=input$comment,pt.size = 2)+
				labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")+guides(colour = guide_legend(title = "Cell type"))+
				scale_color_manual(values = getPalette(colourCount))
		#print(p1)
		#选择颜色
		if(unique(d$Group)=="Drug_resistance"){color="#1b9e77"}else{color="#d95f02"}
		p2 <- DimPlot(d,label = FALSE,group.by = "Group",reduction=input$comment, cols = color, pt.size =1)+
        labs(title="")+xlab("UMAP_1")+ylab("UMAP_2")
		#p2 <- DimPlot(d, label = FALSE,reduction   = input$comment,group.by="Group")
		print(p1+p2)
			}
		 },
		  #输出图片的长宽
		  width = function() input$plotWidth, height = function() input$plotHeight
		)
		# VlnPlot
		output$VlnPlot<-renderPlot({
			c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
			if(input$ChooseSample=="All"){
			VlnPlot(c, features = input$gene, pt.size = 0.2, ncol = 1)
			}else{
			d <- subset(c,sample==input$ChooseSample)
			VlnPlot(d, features = input$gene, pt.size = 0.2, ncol = 1)
			}
			#VlnPlot(pb, features = input$gene, pt.size = 0.2, ncol = 1)
		  },
		  #输出图片的长宽
		  width = function() input$plotWidth, height = function() input$plotHeight
		  )
        #FeaturePlot		  
		output$FeaturePlot<-renderPlot({
			c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
			if(input$ChooseSample=="All"){
			FeaturePlot(c, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
			}else{
			d <- subset(c,sample==input$ChooseSample)
			FeaturePlot(d, features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
			}
			#FeaturePlot(IdentRename(pb,input$Choosecluster,input$clustername), features = input$gene, reduction= input$comment, pt.size = 0.2, ncol = 1)
		  },
		  #输出图片的长宽
		  width = function() input$plotWidth, height = function() input$plotHeight
		  )	  
		#DotPlot	  
	    output$DotPlot<-renderPlot({
			c <- IdentRename(pb,input$Choosecluster,input$Choosecluster)
			gene_names <- strsplit(input$geneinput, "\n")[[1]]
			if(input$ChooseSample=="All"){
			DotPlot(c, features = rev(as.character(unique(gene_names))))+ 
			theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
			ylab("Cell type")+
			xlab("Marker gene")+
			coord_flip()+
			scale_color_viridis()
			}else{
			d <- subset(c,sample==input$ChooseSample)
			DotPlot(d, features = rev(as.character(unique(gene_names))))+ 
			theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
			ylab("Cell type")+
			xlab("Marker gene")+
			coord_flip()+
			scale_color_viridis()
			}
		  },
		  #输出图片的长宽
		  width = function() input$plotWidth, height = function() input$plotHeight
		  )	  
			 		  
	  })
		
	#Metadata 模块
   ### table 1 metadata
  # output$table1 <- DT::renderDataTable({table1},
                                              # options = list(pageLength = 15),
                                              # selection = list(mode='single', selected=NULL)				  
  # )
   
  # Filter data based on selections
  #SC_lung
  output$table1 <- DT::renderDataTable(DT::datatable({
    data <- pb@meta.data[,c(1:6,36,38:41,46,54:56)]
    if (input$sample1 != "All") {
      data <- data[data$sample == input$sample1,]
    }
    if (input$celltype1 != "All") {
      data <- data[data$celltype.main == input$celltype1,]
    }
    if (input$Group1 != "All") {
      data <- data[data$Group == input$Group1,]
    }
    data
	}),options = list(pageLength = 15)#,selection = list(mode='single', selected=NULL)
	
  )
  
    output$table2 <- DT::renderDataTable(DT::datatable({
    data <- pb@assays$RNA@meta.features #LR
    if (input$pathway1 != "All") {
      data <- data[data$pathway_name == input$pathway1,]
    }
    if (input$Group2 != "All") {
      data <- data[data$group == input$Group2,]
    }
    data
	}),options = list(pageLength = 15),selection = list(mode='single', selected=NULL)
	
  )
  #NSCLC
  # output$table3 <- DT::renderDataTable(DT::datatable({
    # data <- nsclc@meta.data#[,c(1:6,36,38:41,46,54:56)]
    # if (input$sample1 != "All") {
      # data <- data[data$sample == input$sample1,]
    # }
    # if (input$celltype1 != "All") {
      # data <- data[data$celltype.main == input$celltype1,]
    # }
    # if (input$Group1 != "All") {
      # data <- data[data$Group == input$Group1,]
    # }
    # data
	# }),options = list(pageLength = 15)#,selection = list(mode='single', selected=NULL)
	
  # )
  
    # output$table4 <- DT::renderDataTable(DT::datatable({
    # data <- nsclc@assays$RNA@meta.features #LR
    # if (input$pathway1 != "All") {
      # data <- data[data$pathway_name == input$pathway1,]
    # }
    # if (input$Group2 != "All") {
      # data <- data[data$group == input$Group2,]
    # }
    # data
	# }),options = list(pageLength = 15),selection = list(mode='single', selected=NULL)
	
  # )
  
  
##### sankey
  # output$sankey <- renderUI({
    # includeHTML("data/IV_immu.html")
  # })

# ##### codeocean result.html
  # output$result <- renderUI({
    # includeHTML("data/result.html")
  # })
  
  # ######添加搜索按钮-----------
  
  # observeEvent(input$cancelButton, {
    # # 在这里添加清空textInput的逻辑
    # shinyjs::runjs('$("#searchInput").val("");')
  # })
  
  # observeEvent(input$searchButton, {
    # # 在这里添加搜索逻辑
    # keyword <- input$searchInput
    # selected_cancer <- input$var
    # # 这里可以执行搜索操作，例如查询数据库或过滤数据
    # # 这里只是一个简单的演示，显示搜索关键字和所选癌症
    # output$searchResult <- renderText({
      # paste("搜索结果：", keyword, "，所选癌症：", selected_cancer)
    # })
  # })
  
  # #设置下载按钮
  # output$Download <- downloadHandler(
    # filename = function() {
      # paste0(input$dataset, ".tsv")
    # },
    # content = function(file) {
      # vroom::vroom_write(data(), file)
    # }
  # )
  
  #download下载数据
   # Reactive value for selected dataset ----
  datasetInput <- reactive({
    switch(input$dataset,
           # "rock" = rock,
           # "pressure" = pressure,
           # "cars" = cars)
		   "LUAD" = pb)#,
		   # "NSCLC" = nsclc)
		   
  })
   # Reactive value for selected datatype ----
  datatypeInput <- reactive({
    switch(input$datatype,
		   "Meta data" = data.frame(datasetInput()@meta.data),
		   # "Expression Matrix" = data.frame(datasetInput()@assays$RNA@counts),
           "Significant Ligand-Receptor Interaction Pairs" = data.frame(datasetInput()@assays$RNA@meta.features)
		   # "Meta data" = data.frame(pb@meta.data),
		   # "Expression Matrix" = data.frame(pb@assays$RNA@counts),
           # "Significant Ligand-Receptor Interaction Pairs" = data.frame(pb@assays$RNA@meta.features)
	)	   
  })

  # Table of selected dataset ----
  #没有表格
  # output$table <- renderTable({
    # datasetInput()
  # })
  
  #有表格展示
  output$table <- DT::renderDataTable({
	df <- datatypeInput()
    action <- DT::dataTableAjax(session, df, outputId = "table")
    DT::datatable(df)#, options = list(pageLength = 15))
  })

  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      #paste(input$dataset, ".csv", sep = "") #只设置输出csv时
	  paste(input$dataset, input$filetype, sep = ".")
    },
    content = function(file) {
	  #write.csv(datasetInput(), file,row.names = FALSE)#只设置输出csv时
	  sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      write.csv(datasetInput(), file, sep=sep,row.names = FALSE)
    }
  )
  
 
}

# Run the application 运行app
shinyApp(ui = ui, server = server)


