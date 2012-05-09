function debug
% print information about Matlab setup
strcat('Userpath: ' , userpath)

strcat('version: ',version)
strcat('release date: ',version('-date'))

which cmu.units
