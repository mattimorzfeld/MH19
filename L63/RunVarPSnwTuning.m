function RunVarPSnwTuning(Nes,SetUps)

inflAll = .9:.05:2;
for jj=1:length(Nes) %[20 50 100 200 500 1000]
    Ne = Nes(jj);
    for kk=1:length(SetUps)%1:4
        SetUp = SetUps(kk);
        for Gap = [2 4 6 8 10 12 14]
            TuneVarPSnw(Ne,inflAll,SetUp,Gap);
        end
    end
end    
