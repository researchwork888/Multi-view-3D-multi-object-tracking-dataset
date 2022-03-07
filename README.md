# Multi-view-3D-People-Tracking-Dataset

CMC Dataset - Multi-view 3D People Tracking Dataset

This 4-camera dataset is published in relation to the following paper: 

@article{ong2020bayesian,
  title={A Bayesian Filter for Multi-view 3D Multi-object Tracking with Occlusion Handling.},
  author={Ong, J and Vo, BT and Vo, BN and Kim, DY and Nordholm, S},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year={2020}
}

Url: https://ieeexplore.ieee.org/document/9242263


Description:

1) CMCX -> Folders containing the respective datasets, i.e.,\
- CMC1 has a maximum of 3 people and virtually no occlusion. \
- CMC2 has a maximum of 10 people with some occlusion. \
- CMC3 has a maximum of 15 people with significant occlusion. \
- CMC4 involves people jumping and falling with a maximum of 3 people. \
- CMC5 involves people jumping and falling with a maximum of 7 people. 

2) MOT format annotations
- Reference: https://motchallenge.net/instructions/
- GT_CMCX_Cam_X.txt -> Bounding-box ground truths for each camera. 
- GT_CMCX_WOLD_CENTROID -> 3D ground truths. 
- DET_CMCX_CamX.txt -> Bounding-box detections for each camera. 



