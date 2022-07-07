# FOGtool #
The golden standard for freezing of gait (FOG) assessment is video annotation by two independent raters. The FOGtool implements the theoretical guideline described by Cockx, Klaver and colleagues[^1], to combine the annotions of two raters based on two parameters: the tolerance and correction parameters. Moreover, it calculates the agreement between the two raters by providing the positive agreement, the negative agreement, and the prevalence index.

In summary, videos are first annotated by two raters as proposed by Gilat in the open-source ELAN software, including characterizing the phenotype (i.e. trembling, shuffling, or akinesia) and trigger (e.g. FOG_Target, FOG_180_R, FOG_Doorway) of each FOG episode [^2]. The annotations are subsequently exported as tabular files which are read in by the FOGtool. The FOG tool will then, based on the chosen tolerance and correction parameter, define the agreed epochs (FOG or no FOG) and the epochs that should be reviewed. Furthermore, FOG episodes that were annotated by both raters, but were characterized by a different phenotype or trigger, are flagged by a ‘check_type’ or ‘check_trigger’, respectively. The outcome is visualized and exported as tabular files which can be reimported into ELAN. Hence, researchers can discuss the remaining to-be-discussed areas (to keep or remove) and the phenotype/trigger of the episode while reviewing the videos. 

## Copyright and Citation ##
FOGtool is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file LICENSE for more details.

Please refer to this paper when using FOGtool [paper is to be included here once accepted].

## Installation ## 
To start using FOGtool, download the latest version [here](https://github.com/helenacockx/FOGtool/releases/latest/) under *assets*.
Download the FOGtool.exe file and follow the download instructions. <br>
When downloading the Source Code.zip, you also have access to example files (in the folder *examplefiles*), and the original functions (combine_FOGannotations.m and agreement_calculator.m in the folder *code*) to start building your own scripts.

A clear user manual can be found [here](https://github.com/helenacockx/FOGtool/blob/main/manual/UserManual.md).

## Feedback ##
Please note that our tool is still new and is open to further improvements. We are happy to receive your feedback and invite you to share your issues or feature enhancements on github or send us an email (helena.cockx@donders.ru.nl or emilieklaver@mst.nl).

[^1]: Cockx H., et al., The Grey Area of Freezing Annotation: a Guideline and Open Source Practical Tool, preprint: https://osf.io/274ch
[^2]: Gilat M, How to annotate freezing of gait from video: A standardized method using open-source software, Journal of Parkinson’s Disease, Accepted for publication on 12th August 2019, DOI: 10.3233/JPD-191700.
