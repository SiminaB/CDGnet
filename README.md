# CDGnet

Repository for code related to the CDGnet tool for prioritizing targeted cancer therapies

The CDGnet tool is currently hosted at http://epiviz.cbcb.umd.edu/shiny/CDGnet/. It is an informatics tool for recommending targeted therapies to individuals with cancer using biological networks. Note that it is a research tool and should not be used for clinical care.
This tool is described in the preprint https://www.biorxiv.org/content/10.1101/605261v1 by Kancherla, Rao, et al.

Please open an issue at https://github.com/SiminaB/CDGnet/issues or email smb310@georgetown.edu if you have any questions or concerns about this software.

This code is distributed under the GPL (>= 3) license. Please note that we are including data from DrugBank and that, while "DrugBank is offered to the public as a freely available resource," commercial use requires a license and
users who use CDGnet with the DrugBank data should also cite the DrugBank papers. More information is available [here](https://www.drugbank.ca/).
**This software comes without warranty of any kind. It is intended for research purposes only and not for clinical care.**

## Repository organization

#### Code of conduct
* _Code of conduct.md_: Required reading for developers and individuals who open issues or make pull requests related to the CDGnet tool.
#### CDGnet documentation
* _CDGnet\__documentation.pdf_: Provides information on using the online CDGnet tool. In order to use the R code directly, please look at the example in the code_notebook directory and the source code in the code directory.
#### Contribution guidelines
* *Contribution guidelines.md*: Document providing details for individuals interested in collaborating, contributing, or asking questions about CDGnet. Please read it before posting issues or making pull requests.
#### Contributors
* *Contributors.md*: List of CDGnet contributors. Individuals who make significant contributions to CDGnet may make a pull request to be added to this list.
#### code directory
* *code/app.R*: Code to create the CDGnet tool at http://epiviz.cbcb.umd.edu/shiny/CDGnet/.
* *code/functions.R*: Functions used in code/app.R. The main wrapper functions are get_cat_1_2 for obtaining targeted therapies from the first 2 categories, as described in the 
preprint (approved either for the alterations in the tumor type or in other tumor types) and get_cat_3_4 for obtaining targeted therapies from categories 3 and 4 (using a network
approach to prioritize additional therapies)
* *code/fda_parser.R*: Code used to obtain the list of FDA-approved drugs. In order to run it, create a "data" subdirectory and download the Products.txt and Submissions.txt files from the zipped folder at https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files.
#### code_notebook directory
* *code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd*: Rmd file that can be used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
#### code_process_KEGG directory
* *code_process_KEGG/get_all_KEGG_hsa_pathways.Rmd* and *.html*: Rmd file that downloads all the KEGG hsa (Homo sapiens) pathways as KGML files (html is output file)
* *code_process_KEGG/parse_KGML.Rmd* and *.html*: Rmd file that parses KGLM files (html is output file)
* *code_process_KEGG/parse_KEGG_info_oncogenes.Rmd* and *.html*: Rmd file that parses oncogene information from KEGG (html is output file)
* *code_process_KEGG/Preprocess_KEGG_objects.Rmd* and *.html*: Rmd file that preprocesses KEGG objects obtained from running files above in order to use standard gene names and disease names (html is output file)
#### data directory (example input tsv and csv files)
* *data/example_input.tsv*: Tab-delimited example file used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
* *data/example_input.csv*: Comma-delimited example file used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
#### objs directory
* *objs/example_CDGnet_database_inputs.RData*: R objects used in the online CDGnet tool as well as the code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd file. In order to be able to run code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd, need additional objects that can be obtained
by running the code in code_process_KEGG


