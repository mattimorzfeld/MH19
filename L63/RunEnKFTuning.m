%% 
clear
close all
clc

inflAll = .9:.05:2;
for Ne = [20 50 100 200 500]
    for SetUp = 5%1:4
        for Gap = [2 4 6 8 10 12 14]
            TuneEnKF(Ne,inflAll,SetUp,Gap);
        end
    end
end
