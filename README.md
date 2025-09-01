p2p-Cortical: A model of cortical electronic prosthetic technologies
====================================================================
This code implements computational model or ‘virtual patient’, based on the neurophysiological architecture of visual cortex, which successfully predicts the perceptual experience of participants across a wide range of previously published human cortical stimulation studies describing the location, size, brightness and spatiotemporal shape of electrically induced percepts in humans. 

Simulations such as the above are likely to be critical for providing realistic estimates of prosthetic vision, thus providing regulatory bodies with guidance into what sort of visual tests are appropriate for evaluating prosthetic performance, and improving current and future technology.

If you use p2p-cortical in a scholarly publication, please cite as:

Fine, I., Boynton, G.M., A virtual patient simulation modeling the neural and perceptual effects of human visual cortical stimulation, from pulse trains to percepts. Scientific Reports, 2024 14, 17400.
[https://www.nature.com/articles/s41598-024-65337-1#citeas](https://doi.org/10.1038/s41598-024-65337-1)

Most of the code in this library was originally written by Ione Fine and Geoffrey Boynten and reviewed, edited, added to and cleaned up by Eirini Schoinas.

Getting started
===============

The bulk of the routines in p2p-cortical are in the Methods file p2p_c.
Example files to get you started are described in the site wiki as well. Additionally, many of the files produce figures from the paper, but can be adapted for other things.

Figures in the paper
===============

To produce the figures in the paper use the following files. You may get slightly different figures due to the inherent randomness in the modeling or also possibly slight tweaks made to the models.

## Figure 1
To be updated. Not in GitHub repo yet.  

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig1_HTML.jpg" width="400">

## Figure 2
Run SimulateMaps.m to recreate these figures.  
Note that some small editions were added in photoshop. 

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig2_HTML.jpg" width="400">

## Figure 3
**3A and 3B:** Run Fit_Freq_and_Pulse_Width.m   
**3C:** Run Winawer_Brightness_Neuron2016.m

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig3_HTML.jpg" width="400">

## Figure 4
**4A, 4B, 4C:** Run Winawer_Size_Neuron2016.m  
**4D:**  Left panel replotted from 
Bosking, W. H., Beauchamp, M. S. & Yoshor, D. Electrical stimulation of visual cortex: Relevance for the development of visual cortical prosthetics. Annu. Rev. Vis. Sci. 3, 141–166. https://doi.org/10.1146/annurev-vision-111815-114525 (2017).  
Right Panel shows corresponding simulations for two eccentricities. Bosking_SizeAmplitude_JNeuro2017.m

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig4_HTML.jpg" width="400">

## Figure 5
**5A:**  Anatomical panel from Neuron, 92/6, J. Winawer and J. Parvizi, linking electrical stimulation of human primary visual cortex, size of affected cortical area, Neuronal Responses, and Subjective Experience, Fig. 1A, Copyright (2016), with permission from Elsevier. 

The white panels show single typical phosphene replotted from from Winawer, J. & Parvizi, J. Linking electrical stimulation of human primary visual cortex, size of affected cortical area, neuronal responses, and subjective experience. Neuron 92, 1213–1219. https://doi.org/10.1016/j.neuron.2016.11.008 (2016). Run Winawer_PhosphenePictures_Neuron2016.m to get the black panels.

**5B:** Run Bosking_SizeEccentricity_JNeuro2017.m

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig5_HTML.jpg" width="400">

## Figure 6
**6A:** Replotted from Beauchamp, M. S. et al. Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell 181, 774-783e775. https://doi.org/10.1016/j.cell.2020.04.033 (2020).  
**6B:** Run BeauchampFig4.m  
**6C:** The upper left panel replots the patient reported phosphene maps of stimulated electrodes (bold circles) and the direction of the temporal sequence of stimulation (arrow) from Beauchamp et al. The lower left panel replots the participant’s actual drawing of the visual percept from Beauchamp et al. Run Beauchamp_DynamicLetters_Cell2020_RF.m and Beauchamp_DynamicLetters_Cell2020_movie.m to get simulated letters in the right panel.

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig6_HTML.jpg" width="400">

## Figure 7
**7A:** Array image and informal observations kindly supplied by P. Troyk and G. Dagnelie. For the rest of the images run Simulate_WFMA.m  
**7B and 7C** Run Phosphene_Sz_vs_Ecc_Simulate. Note that to get the specifc phosphenes you will want to only simulate for one eccentricty/location.

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig7_HTML.jpg" width="400">

## Figure 8
**For 8A, 8B, and 8C:**
Run Array_Sim_call.m which calls Array_Sim_main.m and Array_Sim_movie.m Array_Sim_main.m creates the first two panels. Array_Sim_movie creates the third. The cat image is a frame from CatNotCooperating.m4v or CatNotCooperating.avi.

<img src="https://github.com/ionefine/p2p-cortical/blob/main/wiki_pictures/41598_2024_65337_Fig8_HTML.jpg" width="400">




