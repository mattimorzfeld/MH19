%% 
clear
close all
clc

inflAll = .9:.05:2;
for Ne = [100 200 500 1000 5000]
    for SetUp = 5%1:4
        for Gap = [2 4 6 8 10 12 14]
            TunePFGR(Ne,inflAll,SetUp,Gap);
        end
    end
end
