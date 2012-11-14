portal-rodent-dispersal
=======================

Research on rodent dispersal patterns near Portal AZ. Contains code for reproducing the results of analyses using individual-level tag dat for common rodent species at the Portal Project site. Collaborators on this project include Sarah R. Supp, S. K. Morgan Ernest, and David Koons. Code was written by Sarah R. Supp.

The code and data in this repository all for the analyses and figures to be fully replicated using a subset of the published Portal dataset which includes individual-level rodent data from 2000-2009. Species evaluated include: Peromyscus eremicus (PE), Peromyscus maniculatus (PM), Onychomys torridus (OT), Onychomys leucogaster (OL), Dipodomys merriami (DM), Dipodomys ordii (DO), Chaetodipus baileyi (PB), Chaetodipus penicillatus (PP), Sigmodon hispidus (SH), Sigmodon fulviventer (SF), Neotoma albigula (NAO), and Reithrodontomys montanus (RM).

NOTE: The project and code in this repository is still under development. 

Requirements: R 2.x and the file containing functions specific to this code, movement_fxns.R.

The analyses can be replicated by changing the working directory at the top of the file stake_movements.R to the location on your computer where you have stored the .R and .csv files.

Code should take approximately XX minutes to run start to finish. Figures should output as pdfs to your working directory. 

Data use: Data is provided in this supplement for the purposes of replication. If you wish to use the data for additional research, they should be obtained from the published source (Ecological Archives E090-118-D1; S. K. Morgan Ernest, Thomas J. Valone, and James H. Brown. 2009. Long-term monitoring and experimental manipulation of a Chihuahuan Desert ecosystem near Portal, Arizona, USA. Ecology 90:1708. doi:10.1890/08-1222.1)

Included Files
stake_movment.R script -- cleans up the data, runs the statistical analyses, and outputs figures.
movement_fxns.R script -- holds the relevant functions for executing the stake_movement.R script.
rodent_data. R  -- 
data.csv -- DATA TO BE ADDED

License
This code is available under a BSD 2-Clause License.

Copyright (c) 2012 Weecology. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information
Sarah Supp's email: sarah@weecology.org and sarah.supp@usu.edu

Sarah's website: http://weecology.org/people/sarahsupp/Sarah_Supp/