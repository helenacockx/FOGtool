# FOGtool User Guide #

## What is the FOGtool? ##
The golden standard for freezing of gait (FOG) assessment is video annotation by two raters. The FOGtool implements the theoretical guideline as described by Cockx, Klaver and colleagues, to combine the annotions of these two raters based on two parameters, the tolerance and correction parameters. Moreover, it calculates the agreement between the two raters by providing the positive agreement, the negative agreement, and the prevalence index.

In summary, videos are first annotated by two raters as proposed by Gilat in the open-source ELAN software, including characterizing the phenotype (i.e. trembling, shuffling, or akinesia) and trigger (e.g. FOG_Target, FOG_180_R, FOG_Doorway) of each FOG episode [^1]. The annotations are subsequently exported as tabular files which are read in by the FOGtool. The tool will then, based on the chosen tolerance and correction parameter, define the agreed epochs (FOG or no FOG) and the to-be-discussed areas. Furthermore, FOG episodes that were annotated by both raters, but were characterized by a different phenotype or trigger, are flagged by a ‘check_type’ or ‘check_trigger’, respectively. The outcome is visualized and exported as tabular files which can be reimported into ELAN. Hence, researchers can discuss the remaining to-be-discussed areas (to keep or remove) and the phenotype/trigger of the episode while reviewing the videos. 

## How should I use the FOGtool? ##
This user guide describes the steps how to use the graphical user interface of the FOGtool, starting from the annotation files from ELAN. However, users are also free to build there own analysis scripts by using the individual functions provided in this repository (combine_FOGannotations.m and agreement_calculator.m).<br>
Example files are provided in the folder *examplefiles*.

### 1.  Creating the annotation files in ELAN ###
To start, we highly recommend to annotate the videos by using the format which was described by Gilat [^1]. The ELAN template provided by Gilat can be adapted to the user's preferences, however, the following tiers are required by the FOGtool and should therefore always be included: the Gait_Task, FOG_Trigger and FOG_Type. All extra tiers will be retained by the FOGtool. For further details for annotating the video, we refer to [the manual provided by Gilat](https://morangilat.com/wp-content/uploads/2019/09/FOG_Scoring_User_guide_v1.0.pdf). <br>
Once video annotations are finished, ensure that the ELAN software is set to English, by clicking Options &rarr; Language &rarr; English (see Figure 1). <br>
Then, to export the data, select File &rarr; Export As &rarr; Tab-delimited Text...(see Figure 2). Then select the export settings as indicated by Figure 3, which is equal to the export format as recommended by Gilat. <br> 
**Important remark: The filename should be of the format ‘participantID_raterID’ (e.g. *sub-PD06_annotations-Helena*, *VS24_EKeng*) or ‘participantID_filename_raterID’ (e.g. *sub-PD06_session-01_annotations-Helena*, *VS24_823_EKeng*), in this particular order, as the FOGtool uses the filename separated by the underscores to read unique file identifiers.**

<figure>
<img src="https://github.com/helenacockx/FOG_annotations/blob/main/manual/Images/ELAN_language.png" style="width:70%">
<figcaption align = "left">Figure 1 - ELAN language settings.</figcaption>
 </figure>
 
<figure>
<img src="https://github.com/helenacockx/FOGtool/blob/main/manual/Images/ELAN_export.png" style="width:50%">
<figcaption align = "left">Figure 2 - Select the ELAN export.</figcaption>
 </figure>
 
 <figure>
 <img src="https://github.com/helenacockx/FOG_annotations/blob/main/manual/Images/ELAN_Export_options.PNG" style="width:80%">
<figcaption align = "left">Figure 3 - ELAN export settings.</figcaption>
 </figure>
 
### 2. Running the FOGtool ###
Prior to running the tool, make sure that you have all the files that you want to compare within one run, are organized in two folders: one folder for each rater. <br>
An overview of the functionalities of the FOGtool is given in Figure 4.<br>
Click on button **1a** to select the annotation files of rater 1; then click button **1b** to select the annotation files of rater 2. The selected files should be listed in the preview screen below the buttons. 
Next, choose the values of the tolerance **(2)** and correction **(3)** parameter that should be used to define whether not-overlapping annotation parts should be included as FOG, excluded as FOG, or marked to be further discussed. Not-overlapping parts that are exceding the tolerance parameter (the default is set to 2 seconds) will be marked as 'to discuss'. Not-overlapping parts that are shorter than the tolerance parameter are included/excluded based on the correction parameter (default is 'include'). To learn more about the tolerance and correction parameter, we refer to our paper (Cockx et al. *link to paper is to be included here*).<br>
To determine where the outcome files after running the FOGtool will be exported to, click on button **4** and select the desired folder.<br>
From all the files, the agreement info will be collected in an agreement table which will be used to calculate the overall agreement of a project. Click on button **5** to determine the name and location of this agreement table (default is 'agreement_table.tsv'). Notice that, if you do not run all files at once, you can choose to load a pre-existing agreement table. This will add the new files to the existing table.

Finally, click RUN (button **7**) to calculate the output files and the agreement info of the selected files. This will also displays an image with a graphical representation of the annotations of rater 1, the annotations of rater 2 and the outcome for each file (Figure 5). If required, the figures can be saved by checking the checkmark (functionality 6). To calculate the agreement of the all the files that are saved in the agreement table, click button **8**. 

<figure>
 <img src="https://github.com/helenacockx/FOGtool/blob/main/manual/Images/FOGtool_interface.PNG" style="width:100%">
<figcaption align = "left">Figure 4 - The FOG tool.
 
1.	To select the annotation files of both rater 1 and 2.
2.	To set the tolerance in seconds.
3.	To set the correction parameter to include or exclude.
4.	To select the annotation output folder.
5.	To select the filename and location of the agreement table. 
6.	Check to save the figures which display the annotations of rater 1, the annotations of rater 2, and the outcome of each file.
7.	To run the FOGtool for the selected files.
8.	To calculate the final agreement over all files based on the loaded agreement table. </figcaption>
 </figure>
 
<figure>
 <img src="https://github.com/helenacockx/FOGtool/blob/main/manual/Images/FOGtool_result.PNG" style="width:60%">
<figcaption align = "left">Figure 5 - Example of graphical representation of annotations. </figcaption>
 </figure> 

### 3. Review the combined annotation files in ELAN ###
To evaluate the annotations which were not in full consensus, create a new file in ELAN containing the video file you want to review. 
Add the results from the FOG tool by clicking File &rarr; Import &rarr; CSV / Tab-delimited Text file... and select the .txt file that corresponds to the video. 
Click Open and OK, and the annotations are loaded into ELAN. An example is given in Figure 6. The agreed annotations are indicated by the tiers FOG_agreed_Trigger and FOG_agreed_type. The check_annotation tier indicates FOG episodes that were annotated by both raters, but that contained a different label for the FOG type or FOG trigger. The FOG_disagreed_Trigger and FOG_disagreed_Type tiers indicate the FOG episodes that should be discussed with a third rater, or untill consensus is reached. The output of this discussion can again be exported as described here above, and further processed depending on the study goal.<br>
**Note that neighbouring agreed and disagreed FOG parts might be exported as two separate events, while they belong to the same episode. We therefore recommend the use of % time frozen instead of the number of FOG episodes as primary outcome; or to combine both events into one episode before further processing.**

<figure>
 <img src="https://github.com/helenacockx/FOGtool/blob/main/manual/Images/ELAN_result%20-%20kopie.png" style="width:80%">
<figcaption align = "left"> Figure 6 - The results from the FOG tool in ELAN.</figcaption>
 </figure>
 
[^1]: Gilat M, How to annotate freezing of gait from video: A standardized method using open-source software, Journal of Parkinson’s Disease, Accepted for publication on 12th August 2019, DOI: 10.3233/JPD-191700.
