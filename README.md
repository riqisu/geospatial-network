# geospatial-network
Source code for paper "Data-based reconstruction of complex geospatial networks, nodal positioning and detection of hidden nodes";

This project uses compressive sensing to analyze time series of spatial distributed nodes in network. We use compressive sensing to infer the delay coupling between all nodes, and further reconstruct their spatial location.

Abstract for the paper reads:
Given a complex geospatial network with nodes distributed in a two-dimensional region of physical space, can the locations of the nodes be determined and their connection patterns be uncovered based solely on data? We consider the realistic situation where time series/signals can be collected from a single location. A key challenge is that the signals collected are necessarily time delayed, due to the varying physical distances from the nodes to the data collection centre. To meet this challenge, we develop a compressive-sensing-based approach enabling reconstruction of the full topology of the underlying geospatial network and more importantly, accurate estimate of the time delays. A standard triangularization algorithm can then be employed to find the physical locations of the nodes in the network. We further demonstrate successful detection of a hidden node (or a hidden source or threat), from which no signal can be obtained, through accurate detection of all its neighbouring nodes. As a geospatial network has the feature that a node tends to connect with geophysically nearby nodes, the localized region that contains the hidden node can be identified.

This is a open-access paper and can be downloaded from Royal Society of Science:
http://rsos.royalsocietypublishing.org/content/3/1/150577
