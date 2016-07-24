# Installing necessary dependencies

library(shiny)
library(YavAch)
library(dplyr)
library(RColorBrewer)
library(vcd)
library(stringr)
library(GGally)
library(DT)
library(reshape2)
library(cluster)
library(corrplot)

shinyUI(fluidPage(

##---------------------------------------------------------------------------------------------------------------------

navbarPage(theme = "bootstrap.min4.css",
          title = "Achilles Data Analysis",

#---------------------------------------------------------------------------------------------------------------------

tabPanel("Introduction to the data : RNAi",
                      titlePanel("Introduction to the data : RNAi"),
                      sidebarLayout(
                        sidebarPanel(title = "What is RNAi?", includeMarkdown("Intro_1.Rmd")),
                        mainPanel(
                          tabsetPanel(
                            tabPanel(title = "shRNA mechanism", includeMarkdown("Intro_2.Rmd")),
                            tabPanel(title = "Screens and Limitations", includeMarkdown("Intro_3.Rmd")),
                            tabPanel(title = "Project Achilles", includeMarkdown("Pre_1.Rmd"), plotOutput("cancersPieChart"))
                          )))),

##--------------------------------------------------------------------------------------------------------------------

navbarMenu("shRNA Targetting",

#---------------------------------------------------------------------------------------------------------------------
tabPanel("shRNA Variation Within Individual Cell Lines",
                    titlePanel("shRNA Variation Within Individual Cell Lines"),
                    sidebarLayout(

                      sidebarPanel(
                        DT::dataTableOutput("cancersList"),

                        textInput("geneName_variation.within.cell.line", label = "Gene Name", "MDM4"),

                        submitButton(text = "Apply Changes", icon("refresh")),

                        hr(),
                        includeMarkdown("shT_1.Rmd")

                        ),

                      mainPanel(plotOutput("variation.within.cell.line.plot", height = 900))

                    )),

tabPanel("shRNA Variation Within Cancer Entities",
         titlePanel("shRNA Variation Within Cancer Entities"),
         sidebarLayout(
           sidebarPanel(

             textInput("geneName_variation.within.entity", label = "Gene Name", "MDM4"),

             selectInput("entityName_variation.within.entity", label = h4("Entity"),

                         choices = list("Bone (6)" = "BONE", "Breast (13)" = "BREAST",
                                        "Central Nervous System (35)" = "CENTRAL_NERVOUS_SYSTEM", "Endometrium (2)" = "ENDOMETRIUM", "Haematopoietic and Lymphoid Tissue (30)" = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                                        "Kidney (10)" = "KIDNEY", "Large Intestine (20)" = "LARGE_INTESTINE", "Liver (1)" = "LIVER",
                                        "Lung (21)" = "LUNG", "Oesophagus (10)" = "OESOPHAGUS", "Ovary (29)" = "OVARY",
                                        "Pancreas (17)" = "PANCREAS", "Pleura (2)" = "PLEURA", "Prostate (3)" = "PROSTATE", "Skin (7)" = "SKIN", "Small Intestine (1)" = "SMALL_INTESTINE",
                                        "Soft Tissue (2)" = "SOFT_TISSUE", "Stomach (4)" = "STOMACH",
                                        "Urinary Tract (3)" = "URINARY_TRACT")),

             submitButton(text = "Apply Changes", icon("refresh")),

             hr(),
             includeMarkdown("shT_2a.Rmd"),
             hr(),

             includeMarkdown("shT_2b.Rmd")),

           mainPanel(
             tabsetPanel(
               tabPanel("Option 1", plotOutput("variation.within.entity.plot", height = 800)),
               tabPanel("Option 2", plotOutput("variation.within.entities.plot", height = 900))
                        )))),

tabPanel("Consistency of shRNA Phenotypes",
         titlePanel("Consistency of shRNA Phenotypes"),
         sidebarLayout(
           sidebarPanel(

             textInput("geneName_consistency.of.phenotypes", label = "Gene Name", "MDM4"),
             selectInput("entityName_consistency.of.phenotypes", label = h4("Cancer"),

                         choices = list("Bone (6)" = "BONE", "Breast (13)" = "BREAST",
                                        "Central Nervous System (35)" = "CENTRAL_NERVOUS_SYSTEM", "Endometrium (2)" = "ENDOMETRIUM", "Haematopoietic and Lymphoid Tissue (30)" = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                                        "Kidney (10)" = "KIDNEY", "Large Intestine (20)" = "LARGE_INTESTINE", "Liver (1)" = "LIVER",
                                        "Lung (21)" = "LUNG", "Oesophagus (10)" = "OESOPHAGUS", "Ovary (29)" = "OVARY",
                                        "Pancreas (17)" = "PANCREAS", "Pleura (2)" = "PLEURA", "Prostate (3)" = "PROSTATE", "Skin (7)" = "SKIN", "Small Intestine (1)" = "SMALL_INTESTINE",
                                        "Soft Tissue (2)" = "SOFT_TISSUE", "Stomach (4)" = "STOMACH",
                                        "Urinary Tract (3)" = "URINARY_TRACT")),

             submitButton(text = "Apply Changes", icon("refresh")),

             includeMarkdown("shT_3.Rmd")
             ),

             mainPanel(
             tabsetPanel(
               tabPanel("All" , plotOutput("consistency.of.phenotypes.all.plot", height = 800)),
               tabPanel("By Chosen Entity" , plotOutput("consistency.of.phenotypes.plot", height = 800)))
           )))
),

##--------------------------------------------------------------------------------------------------------------------

navbarMenu("p53",

#---------------------------------------------------------------------------------------------------------------------
tabPanel("p53 status and Chosen Gene",
                    titlePanel("p53 status and Chosen Gene"),
                    sidebarLayout(
                      sidebarPanel(

                        textInput("geneName_p53.dependency.by.entity", label = "Gene", value = "MDM4"),
                        selectInput("entityName_p53.dependency.by.entity", label = h4("Cancer"),
                                    choices = list("Bone (6)" = "BONE", "Breast (13)" = "BREAST",
                                                   "Central Nervous System (35)" = "CENTRAL_NERVOUS_SYSTEM", "Endometrium (2)" = "ENDOMETRIUM", "Haematopoietic and Lymphoid Tissue (30)" = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                                                   "Kidney (10)" = "KIDNEY", "Large Intestine (20)" = "LARGE_INTESTINE", "Liver (1)" = "LIVER",
                                                   "Lung (21)" = "LUNG", "Oesophagus (10)" = "OESOPHAGUS", "Ovary (29)" = "OVARY",
                                                   "Pancreas (17)" = "PANCREAS", "Pleura (2)" = "PLEURA", "Prostate (3)" = "PROSTATE", "Skin (7)" = "SKIN", "Small Intestine (1)" = "SMALL_INTESTINE",
                                                   "Soft Tissue (2)" = "SOFT_TISSUE", "Stomach (4)" = "STOMACH",
                                                   "Urinary Tract (3)" = "URINARY_TRACT"),
                                    selected = "SKIN"),

                        selectInput("status_p53.dependency.by.entity", label = h4("p53 status"),
                                    choices = list("All" = "all", "Mutant" = "MUT", "Wild Type" = "WT", "Unknown" = "NA")),

                        submitButton(text = "Apply Changes", icon("refresh")),
                        hr(),
                        includeMarkdown("p53_1.Rmd")),

                      mainPanel(
                        tabsetPanel(
                          tabPanel("p53 dependent genes by Entity", plotOutput("p53.dependency.by.entity.plot", height = 900)),
                          tabPanel("p53 dependent genes Overall", plotOutput("p53.dependency.by.entity.all.plot", height = 900))
                        )))),

tabPanel("GENE-E Analysis of p53",
         titlePanel("GENE-E Analysis of p53"),
         sidebarLayout(
           sidebarPanel(
             textInput("geneName_GENEE", label = "Gene", value = "MDM4"),
             selectInput("data_GENEE", label = h4("Dataset"),
                         choices = list("Achilles v2.0" = TRUE, "Achilles v2.4" = FALSE)),

             submitButton(text = "Apply Changes", icon("refresh")),
             hr(),
             includeMarkdown("p53_2.Rmd")
           ),

           mainPanel(
             tabsetPanel(
               tabPanel("All Cancers", plotOutput("riger.gene.comparison.all.plot", height = 800)),
               tabPanel("Entity", plotOutput("riger.gene.comparison.entities.plot", height = 800)),
               tabPanel("Achilles v2.0 vs Achilles v2.4.3", plotOutput("riger.gene.comparison.datasets.all.plot", height = 800)),
               tabPanel("Achilles v2.0 vs Achilles v2.4.3 (Entities)", plotOutput("riger.gene.comparison.datasets.entities.plot", height = 800))
             ))))
)
)

))

