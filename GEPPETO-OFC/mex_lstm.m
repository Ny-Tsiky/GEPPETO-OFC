CODEXT = '.cpp';
OBJEXT = '.o';  
LIBEXT = '.a';  

HOME = '/Users/tsikyrakotomalala/Documents/these/GEPPETO-OFC'; % Needs to be modified 
OPTHOME = '.';
INCMYL = ['-I',HOME,'/Includes'];
INCKAL = ['-I',OPTHOME,'/Kalman'];
INCOPT = ['-I',OPTHOME,'/Optimfuns'];
INCLSTM = ['-I','./codegen/codegen/lib/my_predict'];
INCNR = ['-I',HOME,'/Includes/nrd'];
LIBMYL = [HOME,'/Libs/math_mac', LIBEXT];
LIBNR = [HOME,'/Libs/nrd',LIBEXT];
CXXFLAGS = '-g';
res = mex('-c', CXXFLAGS, [OPTHOME, '/Kalman/kalman8c.cpp'], INCMYL, INCNR);
res = res + mex('-c',CXXFLAGS, [OPTHOME, '/Optimfuns/plantkalman_discrete.cpp'], INCMYL, INCNR, INCKAL);
res = res + mex('-c',CXXFLAGS, [OPTHOME, '/Optimfuns/simnoisy_discrete8.cpp'], INCMYL, INCNR, INCKAL); % noisy
res = res + mex('-c', CXXFLAGS, [OPTHOME, '/Optimfuns/dop0viafuns.cpp'], INCMYL,INCNR,INCOPT,INCKAL); %
res = res + mex('-c',CXXFLAGS, 'tongue2form.cpp', INCMYL, INCNR);
res = res + mex('-c',CXXFLAGS, 'tongue_geom_audio_tactile_autoenc.cpp', INCMYL, INCNR);

res = res + mex('-c',CXXFLAGS, 'autoencoder.cpp', INCMYL, INCNR); % 
res = res + mex('-c',CXXFLAGS, 'goal_geometry.cpp', INCMYL, INCNR); % 

% compile once
res = res + mex('-c',CXXFLAGS, [OPTHOME, '/codegen/codegen/lib/my_predict/my_lstmForward_c11.cpp'], INCMYL, INCNR, INCLSTM); % Mex wrapper to tongue2form.cpp
res = res + mex('-c',CXXFLAGS, 'tonguefeedback_3modalities.cpp', INCMYL, INCNR, INCKAL);
res = res + mex('-c',CXXFLAGS, 'tonguefuns_nosmoothing.cpp', INCMYL,INCKAL,INCNR,INCOPT,INCLSTM);

mex('-cxx','-output','lstmmultisens_simnoisy_nosmoothing',CXXFLAGS,['simnoisy_discrete8',OBJEXT],['tonguefuns_nosmoothing',OBJEXT],...
    ['tonguefeedback_3modalities',OBJEXT], ['goal_geometry',OBJEXT],  ...
    ['my_lstmForward_c11',OBJEXT],...
    ['tongue_geom_audio_tactile_autoenc',OBJEXT],['autoencoder',OBJEXT],['plantkalman_discrete',OBJEXT],['kalman8c',OBJEXT], LIBMYL, LIBNR)

mex('-cxx','-output','lstmmultisens_via_dop0funs_nosmoothing',CXXFLAGS,['dop0viafuns',OBJEXT],['tonguefuns_nosmoothing',OBJEXT],...
    ['tonguefeedback_3modalities',OBJEXT], ['goal_geometry',OBJEXT],  ...
    ['my_lstmForward_c11',OBJEXT],...
    ['tongue_geom_audio_tactile_autoenc',OBJEXT],['autoencoder',OBJEXT], LIBMYL, LIBNR, INCLSTM) %,['plantkalman_discrete',OBJEXT],['kalman7',OBJEXT]

