+++
# Date this page was created.
date = 2016-04-27T00:00:00

# Project title.
title = "IBI VizEdit"

# Project summary to display on homepage.
summary = "A brand-new, open-source, RShiny application for processing and editing heart rate data"

# Optional image to display on homepage (relative to `static/img/` folder).
image_preview = "IBI-thumb.jpg"

# Tags: can be used for filtering projects.
# Example: `tags = ["machine-learning", "deep-learning"]`
tags = ["Methods", "IBI-VizEdit"]

# Optional external URL for project (replaces project detail page).
#external_link = ""

# Does the project detail page use math formatting?
math = false

# Optional featured image (relative to `static/img/` folder).
[header]
image = "IBI-Screen.PNG"
caption = "IBI VizEdit's Primary Editing Interface"

+++

IBI VizEdit is a program built using RShiny. It is designed to assist in the manual editing of inter-beat interval files that are derived from photoplethysmogram (PPG) recordings. Unlike the electrocardiogram signal (EKG or ECG), PPG signals are characterized by a slow-moving waveform, which presents a different set of challenges when the true signal becomes corrupted by motion artefacts and other sources of noise.

Though increasingly popular due to their ease of use, most heart rate editing software that exists to date was designed and optimized for the detection and editing of inter-beat interval files derived from ECG signals. IBI VizEdit provides a new suite of tools for researchers who find themselves working with messy PPG files.

Please note that IBI VizEdit is beta software. It has not been fully tested, and there are likely numerous bugs and opportunities to optimize code and performance. Any and all feedback is welcome. 

As of right now, IBI VizEdit is only supported for use on Windows 7/8/10 and Linux (Ubuntu 16.04 in particular).

Please cite as:

Barstead, M. G. (2018). IBI VizEdit v.1.2-beta: An RShiny Application [Computer software]. University of Maryland. doi: 10.5281/zenodo.1209474