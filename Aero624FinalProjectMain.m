% Aero 624 Final Project
% The overall goal of this project is too code a simulation that relicates
% the J2 effect that earth has on a sattelite's orbital elements in MATLAB, 
% and then produce the same result using the Unreal Engine. UE4 is the
% platform that our research group, SpaceCRAFT, advised by Dr. Chamitoff, 
% uses to code our simulations.

% I am starting with MATLAB because I have become very familiar with the
% syntax and nuances of using the built in numerical solvers like ode45.
% All of the simulations within the unreal engine are coded in C++ and
% don't have such luxuries, so a large task of the project will be
% learning how to solve differential equations in C++.

% With that said, there are many possible applications to including this
% functionality to SpaceCRAFT and it's why I've chosen this particular
% task. If a user wants to simulate how a particular orbit configuration
% will affect their satellite, they can simply "plug in" the functionality
% of the code and the effect that J2 has can be analyzed in the scope of
% their mission.

% The format of the code will be set up to where the user can specify the
% mission and satellite parameters, and then visualize what the orbit of
% the satellite will look like over time. The initial tests will use orbits
% with a high enough altitude that drag has negligible effects.

% Future iterations of the code will include atmospheric drag effects to
% increase the fidelity of the simulation. I also plan on incorporating the
% third body problem so that users can begin to plan interplanetary
% trajectories. This is a very crucial aspect to the what the
% functionalities should be like iin SpaceCRAFT and why I think it's a good
% option to choose this as my final Project. Do note, that what I turn in
% will not be the final product, work will continue over the summer and
% into the fall with hopes that a robust system can be included in a
% software release at the end of the Fall semester or in the spring.
