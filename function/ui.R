ui <- bs4Dash::dashboardPage(
  ####1.big title####
  title = "Lipid Annotation (beta)",
  bs4Dash::dashboardHeader(title = "Lipid Annotation (NRP V1.0)"),
  ####2.side menu####
  bs4Dash::dashboardSidebar(
    bs4Dash::sidebarMenu(
      ####2.1 part1 sidemenu####
      bs4Dash::menuItem(text = "Introduction",
                        tabName = "part1"),
      ####2.2 part2 sidemenu####
      bs4Dash::menuItem(text = "Peaks processing",
                        tabName = "part2"),
      ####2.3 part3 sidemenu####
      bs4Dash::menuItem(text = "NPG-MS2 match",
                        tabName = "part3"),
      ####2.4 part4 sidemenu####
      bs4Dash::menuItem(text = "RT based PHS",
                        tabName = "part4"),
      ####2.5 part5 sidemenu####
      bs4Dash::menuItem(text = "PG-MS2 match",
                        tabName = "part5"),
      tabItem(
        div(style = 'text-align:center;color:#3c8dbc;font-weight:bold;font-weight: 700;', h4('Console')),
        tabName = "output_tab",
        div(style = "border: 3px dashed #3c8dbc;height:300px;width:220px;",
            verbatimTextOutput("console"))
      )
    )
  ),
  ####3.body####
  bs4Dash::dashboardBody(
    ####3.1 part1 body####
    bs4Dash::tabItem(
      tabName = "part1",
      bs4Dash::blockQuote(
        h1("Introduction", style = 'font-size:36px;'),
        color = 'lightblue',
        style = "height:60px;"
      ),
      tags$p(
        style = "text-align:justify;padding:20px;
             font-family:'Times New Roman',sans-serif;",
        "This method is based on information such as mass spectrometry
             secondary spectrum, relative retention time, ion collision area,
             double bond position, etc., to construct a new non-priori greedy
             secondary matching method ",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>(NPG-MS2 match)</span>"
        ),
        ", combined with relative
             retention time as prior information, Perform a heuristic search",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>(RT based PHS I/II)</span>"
        ),
        ".",
        tags$br(),
        tags$br(),
        "With the help of the relative position information of the retention
             time of different total saturation within the same lipid class,
             the migration association analysis is completed, and the
             high-quality lipid mass spectrum is annotated by combining the
             above three parts. For the low-mass mass spectrum, combined with
             the relative position of its retention time, with the help of
             information enhancement of the related high-quality spectrum,
             the greedy secondary matching ",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>(PG-MS2)</span>"
        ),
        " with prior information is
             completed.",
        tags$br(),
        tags$br(),
        '',
        tags$br()
      ),
    ),
    ####3.2 part2 body####
    bs4Dash::tabItem(
      tabName = "part2",
      bs4Dash::blockQuote(
        h1("Peaks processing", style = 'font-size:36px;'),
        color = 'lightblue',
        style = "height:60px;"
      ),
      tags$p(
        style = "text-align:justify;padding:20px;
             font-family:'Times New Roman',sans-serif;",
        "This part is the preprocessing of mass spectrum peaks, including
             using XCMS and MSnbase packages for centralized processing, reading
             peaks, and peaks processing.",
        tags$br(),
        tags$br(),
        "First, you need to upload the 'mzML' format file of the mass
             spectrum peak, and second, you need to select the ion mode
             (ESImode), where 'P' means positive ion mode and 'N' means
             negative ion mode."
      ),
      div(
        style = "border: 0px;padding:20px;font-family:'Times New Roman',sans-serif;",
        selectInput(
          inputId = "pn_input",
          label = "Choose one ESImode",
          choices = c("P", "N")
        ),
        fileInput(
          inputId = "file",
          label = "Upload the mzML file",
          width = "100%",
          accept = c('.mzML')
        ),
        tags$p(
          style = "text-align:justify;
             font-family:'Times New Roman',sans-serif;",
          "Second, after completing the above operations, wait for the software to
             run until the prompts that ",
          HTML(
            "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  choose peaks processing or not.</span>"
          ),
          "During this process, You need to determine whether the peak needs
             to be processed, and if so, determine the parameters and input them
             in sequence, then press the button. The centralized results will
             be saved as a new 'mzML'
             file, and the sorted results will be saved in a new 'rda' file."
        ),
        fluidRow(
          style = "border: 0px;padding:0px;font-family:'Times New Roman',sans-serif;font-weight: 700;",
          column(
            width = 6,
            "Mass peak centralization (note parameters)",
            selectInput(
              inputId = "smooth_method",
              label = "Choose one smooth method",
              choices = c("MovingAverage", "SavitzkyGolay")
            ),
            numericInput(inputId = "halfWindowSize", "Input half window size", value =
                           0.01),
            selectInput(
              inputId = "pickPeaks_method",
              label = "Choose one pick Peaks method",
              choices = c("kNeighbors", "descendPeak")
            ),
            selectInput(
              inputId = "pickfeature_method",
              label = "Choose one pick Feature method",
              choices = c("centWave", "matchedFilter", "obiwarp", "mzMatch")
            ),
            numericInput(
              inputId = "mzerror",
              "Input Mass-to-charge ratio error (ppm)",
              value = 0.01
            ),
            tags$p('Click the button after check parameters'),
            tags$style(
              HTML(
                "
                   .my-button {
                   background-color: #3c8dbc !important;
                   color: white !important;
                   border-color: transparent !important;
                   box-shadow: 3px 3px 12px rgba(0, 0, 0, 0.5) !important;
                   }
                   "
              )
            ),
            actionButton("start2", "Start Peak processing and picking", class = "btn btn-primary my-button")
          ),
          column(
            width = 6,
            "Mass spectrum peak decentralization",
            tags$p(tags$br()),
            actionButton("start", "Start Peak picking", class = "btn btn-primary my-button")
          )
        )
      ),
    ),
    ####3.3 part3 body####
    bs4Dash::tabItem(
      tabName = "part3",
      bs4Dash::blockQuote(
        h1("NPG-MS2 match", style = 'font-size:36px;'),
        color = 'lightblue',
        style = "height:60px;"
      ),
      tags$p(
        style = "text-align:justify;padding:20px;
             font-family:'Times New Roman',sans-serif;",
        "This part is the MS2 matching algorithm of greedy search without
             prior information.",
        tags$br(),
        tags$br(),
        "First, you need to determine the tolerance of the MS charge-to-mass
             ratio error, then loading the MS2 library,
             wait for the output message ",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  start MS2 result processing.</span>"
        )
      ),
      div(
        style = "border: 0px;padding:20px;font-family:'Times New Roman',sans-serif;",
        numericInput(
          inputId = "ppmerror",
          "Input Mass-to-charge ratio tolerance (ppm)",
          value = 0.01
        ),
        fileInput(
          inputId = "library.file",
          label = "Choose the MS2 library",
          width = "100%",
          accept = c('.rda')
        ),
        tags$p(
          style = "text-align:justify;
             font-family:'Times New Roman',sans-serif;",
          "Second, you need to input matching score threshold, and upload the NPG_MS2
             result. Waiting until the prompt ",
          HTML(
            "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  start RT based PHS.</span>"
          )
        ),
        numericInput(
          inputId = "score",
          "Input matching score threshold (advised value: 0.6)",
          value = 0.01
        ),
        fileInput(
          inputId = "NPG_MS2.file",
          label = "Choose the NPG_MS2.csv file",
          width = "100%",
          accept = c('.csv')
        )
      )
    ),
    ####3.4 part4 body####
    bs4Dash::tabItem(
      tabName = "part4",
      bs4Dash::blockQuote(
        h1("RT based PHS", style = 'font-size:36px;'),
        color = 'lightblue',
        style = "height:60px;"
      ),
      tags$p(
        style = "text-align:justify;padding:20px;
             font-family:'Times New Roman',sans-serif;",
        "In this part, the retention time will be used as prior information
             and combined with the secondary spectrum to further judge and test.",
        tags$br(),
        tags$br(),
        "First, you need to upload the processed file, and wait until the prompt",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  start PSH II.</span>"
        )
      ),
      div(
        style = "border: 0px;padding:20px;font-family:'Times New Roman',sans-serif;",
        fileInput(
          inputId = "processed.file",
          label = "Choose NPG_MS2_processed.csv",
          width = "100%",
          accept = c('.csv')
        ),
        tags$p(
          style = "text-align:justify;
             font-family:'Times New Roman',sans-serif;",
          "Second, you need to upload the PHS I result file, and wait until the prompt",
          HTML(
            "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  start PG-MS2 match.</span>"
          )
        ),
        fileInput(
          inputId = "PHS_I.file",
          label = "Choose PHS_I.csv",
          width = "100%",
          accept = c('.csv')
        )
      )
    ),
    ####3.5 part5 body####
    bs4Dash::tabItem(
      tabName = "part5",
      bs4Dash::blockQuote(
        h1("PG-MS2 match", style = 'font-size:36px;'),
        color = 'lightblue',
        style = "height:60px;"
      ),
      tags$p(
        style = "text-align:justify;padding:20px;
             font-family:'Times New Roman',sans-serif;",
        "This part performs spectrum enhancement for the points whose
             spectrum matches are lower than the set threshold.",
        tags$br(),
        tags$br(),
        "First, you need to upload the peak processed .rda file, and wait until the prompt",
        HTML(
          "<span style='color:#3c8dbc;font-weight:bold;'>Now you can
                  start PSH II.</span>"
        )
      ),
      div(
        style = "border: 0px;padding:20px;font-family:'Times New Roman',sans-serif;",
        fileInput(
          inputId = "peak.rda.data",
          label = "Choose peak processing.rda",
          width = "100%",
          accept = c('.rda')
        ),
        fileInput(
          inputId = "npgms2.data",
          label = "Choose NPG-MS2 result.csv",
          width = "100%",
          accept = c('.csv')
        ),
        fileInput(
          inputId = "part1.data",
          label = "Choose part1 result.csv",
          width = "100%",
          accept = c('.csv')
        ),
        numericInput(
          inputId = "strength.ratio",
          "Input augmentation Index (advised value: 0.2 ~ 0.6)",
          value = 0.01
        ),
        numericInput(
          inputId = "topn",
          "Input maximum candidate",
          value = 1
        ),
        tags$p(style = "border: 0px;font-family:'Times New Roman',sans-serif;font-weight: 700;",
               'Click the button after check parameters'),
        actionButton("start.part2", "Start PG-MS2 match", class = "btn btn-primary my-button")
      ),
      tags$p(
        style = "text-align:justify;
             font-family:'Times New Roman',sans-serif;",
        "Second, upload the PG-MS2 file, to get final output."
      ),
      div(
        style = "border: 0px;padding:20px;font-family:'Times New Roman',sans-serif;",
        fileInput(
          inputId = "pg.data",
          label = "Choose PG-MS2 file",
          width = "100%",
          accept = c('.csv')
        )
      ),
      
      
      
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      ####################################################################################
      ####################################################################################
      
      tags$p(
        style = "text-align:justify;padding:20px;font-weight: 700;font-size:20px;
             font-family:'Times New Roman',sans-serif;",
        HTML(
          "<span style='text-decoration: underline; text-shadow: 1px 1px 1px grey;
                    text-align:justify;font-weight: 700;font-size:25px;
                    font-family:'Times New Roman',sans-serif;'>Citation:</span>"
        ),
        tags$br(),
        "[1] F.Tan, COVID-19, Nature.",
        tags$br(),
        tags$br(),
        HTML(
          "<span style='text-decoration: underline; text-shadow: 1px 1px 1px grey;
                    text-align:justify;font-weight: 700;font-size:25px;
                    font-family:'Times New Roman',sans-serif;'>Acknowledgement:</span>"
        ),
        tags$br(),
        'This project was mainly developed by Hao Xu and Liangsheng Chen,
             supporting by Prof. S.Lin, Prof. J.Zen, Prof. Z.Yang, Z.Hua.,
             and the annotation result was manual checked by Xiaoyun Huang.',
        tags$br(),
        tags$br(),
        HTML(
          "<span style='text-decoration: underline; text-shadow: 1px 1px 1px grey;
                    text-align:justify;font-weight: 700;font-size:25px;
                    font-family:'Times New Roman',sans-serif;'>Connection:</span>"
        ),
        tags$br(),
        "E-mail: xuhaostatistics@163.com; xh852152344@gmail.com",
        tags$br(),
        "Address: National Institute for Data Science in Health and Medicine,
             Xiamen University, Xiamen, Fujian Province, China."
      )
      
      
      
    )
  )
)
