# CDGnet
Repository for code related to the CDGnet tool

The CDGnet tool is currently hosted at http://epiviz.cbcb.umd.edu/shiny/CDGnet/. It is an informatics tool for recommending targeted therapies to individuals with cancer using biological networks. Note that it is a research tool and should not be used for clinical care.
This tool is described in the preprint https://www.biorxiv.org/content/10.1101/605261v1 by Kancherla, Rao, et al.

Please open an issue at https://github.com/SiminaB/CDGnet/issues or email smb310@georgetown.edu if you have any questions or concerns about this software.

This code is distributed under the GPL (>= 3) license.

## Repository organization

* Code of conduct.md: Code of conduct for developers and individuals who comment or open issues related to the CDGnet tool.
* code/app.R: Code to create the CDGnet tool at http://epiviz.cbcb.umd.edu/shiny/CDGnet/.
* code/functions.R: Functions used in code/app.R. The main wrapper functions are get_cat_1_2 for obtaining targeted therapies from the first 2 categories, as described in the 
preprint (approved either for the alterations in the tumor type or in other tumor types) and get_cat_3_4 for obtaining targeted therapies from categories 3 and 4 (using a network
approach to prioritize additional therapies)
* code/fda_parser.R: Code used to obtain the list of FDA-approved drugs. In order to run it, create a "data" subdirectory and download the Products.txt and Submissions.txt files from the zipped folder at https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files.
* code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd: Rmd file that can be used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
* data/example_input.tsv: Tab-delimited example file used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
* data/example_input.csv: Comma-delimited example file used to generate the results presented in the preprint for a putative patient with ER+ breast cancer and FGFR1 overexpression.
* data/example_CDGnet_database_inputs.RData: R objects used in the online CDGnet tool as well as the code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd file. In order to be able to run code_notebook/CDGnet_targeted_therapy_breast_ER_FGFR1.Rmd, need additional objects that can be obtained
by running the code in code_process_KEGG_DrugBank

