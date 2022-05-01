# User guide on how to use the FOG tool #
This user guide describes the steps on how to process annotation files of 2 independent raters, rated by using the ELAN Software, with the FOG Tool. 
## 1. Table of content ##
This one neccessary?
## 2.  Creating the annotation files in ELAN ##
To start annotating the data it is highly recommend to Use the format as described by Gilat [^1]. 
The ELAN template provided by Gilat can be used and altered to suite the to be analysed study, however the FOG Tool requires the following tiers for its processing and should therefore always be included; the Gait_Task, FOG_Trigger and FOG_Type tier. All extra tiers will be retained by the FOG tool. For further details for annotating the data, we refer to [the manual provided by Gilat](https://morangilat.com/wp-content/uploads/2019/09/FOG_Scoring_User_guide_v1.0.pdf). 
Once annotating the data is finished, ensure the ELAN software is set to English, by clicking Options, Language, English, as seen in Figure 1. 
Then, to export the data, select File, Export as and select CSV / Tab-delimited Text file. Then select the export settings as seen in Figure 2, which is equal to the export format as recommended by Gilat. The filename should be of the format ‘filename_raterID’ or ‘subjectname_filename_raterID’, as the FOG tool uses the filename separated by the underscores to read unique file identifiers. 

<figure>
<img src="https://github.com/helenacockx/FOG_annotations/blob/main/Images/ELAN_language.png" style="width:45%">
<figcaption align = "left">Figure 1 - ELAN language settings.</figcaption>
 </figure>
 
 <figure>
 <img src="https://github.com/helenacockx/FOG_annotations/blob/main/Images/ELAN_Export_options.PNG" style="width:45%">
<figcaption align = "left">Figure 2 - ELAN export settings.</figcaption>
 </figure>
 
## 3. Evaluating the annotation files ##
Prior to the analysis, make sure the files you want to compare within one run are in one folder per rater on your computer. Open the FOG tool, an overview of the functionalities is given in Figure 3. 
Click on button 1a to select the annotation files of rater 1, then click button 1b to select the annotation files of rater 2. The selected files should be listed in the preview screen below the buttons. Then, set the tolerance in seconds and correction parameters, as indicated by 3 and 4. The standard tolerance is set to two seconds and the standard correction is to include the rated freezing of gait (FOG)  in the FOG episode.  To learn more about the tolerance and correction, please refer to Cockx et al. **(en dan hier de link naar de paper zelf)**.
To determine where the annotation files after the comparison will be exported to, click on button 4 and select the desired pathway. From all the files an agreement table will be created, click on button 5 to determine the name and location of this agreement table. You notice you can also choose to load a pre-existing agreement table, in order to facilitate that not all data has to be analysed at once. 
To obtain the annotation files containing the ratings of two raters and to calculate the agreement table, click button 7. This also displays an image of the agreement per annotation. If required, the figures can be saved by setting the checkmark indicated by 6. To calculate the agreement of the all the files displayed in the agreement table, click button 8. 

<figure>
 <img src="https://github.com/helenacockx/FOG_annotations/blob/main/Images/FOGtool_nummers_Klad.png" style="width:70%">
<figcaption align = "left">Figure 3 - The FOG tool.
 
1.	To select the annotation files of both rater 1 and 2
2.	To set the tolerance in seconds
3.	To set the correction parameter to include or exclude
4.	To select the annotation output folder
5.	To select the filename and location of the agreement table 
6.	Check to save the figures which display the agreement between raters.
7.	To create the agreement table
8.	To calculate the final agreement over all files. </figcaption>
 </figure>

## 4. Review the combined annotation files in ELAN ##
To evaluated the annotations which were not in full consensus, create a new annotation in ELAN containing the video file you want to review. 
Add the results from the FOG tool by clicking File, Import,  CSV / Tab-delimited Text file and select the .txt file that corresponds to the video. 
Click Open and OK, and the annotations are loaded into ELAN. An example is given in Figure 4. The check_annotation tier indicates which FOG episodes were in consensus, however no consensus was reached on the FOG type or FOG trigger. The FOG_disagreed_Trigger and FOG_disagreed_Type indicate which freezing episodes no consensus was reached. 

<figure>
 <img src="https://github.com/helenacockx/FOG_annotations/blob/main/Images/ELAN_result.png" style="width:70%">
<figcaption align = "left"> Figure 4 - The results from the FOG tool in ELAN.</figcaption>
 </figure>
 
[^1]: Gilat M, How to annotate freezing of gait from video: A standardized method using open-source software, Journal of Parkinson’s Disease, Accepted for publication on 12th August 2019, DOI: 10.3233/JPD-191700.
