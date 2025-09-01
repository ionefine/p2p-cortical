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
Example files to get you started are described in the site wiki.

Figures in the paper
===============

To produce the figures in the paper use the following files:

## Figure 1
To be updated. Not in GitHub repo yet. 

## Figure 2
Run SimulateMaps.m to recreate these figures.  
Note that some small editions were added in photoshop. 

## Figure 3
**3A and 3B:** Run Fit_Freq_and_Pulse_Width.m   
**3C:** Run Winawer_Brightness_Neuron2016.m

## Figure 4
**4A, 4B, 4C:** Run Winawer_Size_Neuron2016.m
**4D:**  Left panel replotted from 
Bosking, W. H., Beauchamp, M. S. & Yoshor, D. Electrical stimulation of visual cortex: Relevance for the development of visual cortical prosthetics. Annu. Rev. Vis. Sci. 3, 141–166. https://doi.org/10.1146/annurev-vision-111815-114525 (2017).  
Right Panel shows corresponding simulations for two eccentricities. Run Phosphene_Sz_vs_Ecc_Simulate.m


## Figure 5
**5A:**  Anatomical panel replotted from Neuron, 92/6, J. Winawer and J. Parvizi, linking electrical stimulation of human primary visual cortex, size of affected cortical area, Neuronal Responses, and Subjective Experience, Fig. 1A, Copyright (2016), with permission from Elsevier. 

The white panels show single typical phosphene replotted from from Winawer, J. & Parvizi, J. Linking electrical stimulation of human primary visual cortex, size of affected cortical area, neuronal
responses, and subjective experience. Neuron 92, 1213–1219. https://doi.org/10.1016/j.neuron.2016.11.008 (2016).

**5B:** Run Bosking

