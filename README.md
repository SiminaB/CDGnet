# CDGnet
Repository for code related to the CDGnet tool

The CDGnet tool is currently hosted at http://epiviz.cbcb.umd.edu/shiny/CDGnet/. It is an informatics tool for recommending targeted therapies to individuals with cancer using biological networks. Note that it is a research tool and should not be used for clinical care.
This tools is described in the preprint https://www.biorxiv.org/content/10.1101/605261v1 by Kancherla, Rao, et al.

Please open an issue at https://github.com/SiminaB/CDGnet/issues or email smb310@georgetown.edu if you have any questions or concerns about this software.

This code is distributed under the GPL (>= 3) license.

## Repository organization

* Code of conduct.md: Code of conduct for developers and individuals who comment or open issues related to the CDGnet tool.
* code/app.R: Code to create the CDGnet tool at http://epiviz.cbcb.umd.edu/shiny/CDGnet/.
* code/functions.R: Functions used in code/app.R. The main wrapper functions are get_cat_1_2 for obtaining targeted therapies from the first 2 categories, as described in the 
preprint (approved either for the alterations in the tumor type or in other tumor types) and get_cat_3_4 for obtaining targeted therapies from categories 3 and 4 (using a network
approach to prioritize additional therapies)
* code/fda_parser.R: Code used to obtain the list of FDA-approved drugs. In order to run it, create a "data" subdirectory and download the Products.txt and Submissions.txt files from the zipped folder at https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files.
