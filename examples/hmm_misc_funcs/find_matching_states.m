function [best_matches, scores, best_scores] = find_matching_states(Gamma_true,Gamma_in)

%  [best_matches, scores, best_scores]=find_matching_states(rand(100,3),rand(100,3)); 

for jj=1:size(Gamma_true,2)
    
    for kk=1:size(Gamma_in,2)
        cc=corrcoef(Gamma_in(:,kk), Gamma_true(:,jj));
        ccs(jj,kk)=cc(1,2);
    end

end

[a b] = max3d(ccs);
best_matches=zeros(size(Gamma_true,2),1);

scores=ccs;

for ii=1:size(Gamma_true,2)
    [jj kk] = max3d(ccs);
    best_matches(jj)=kk;
    ccs(:,kk)=nan;
    ccs(jj,:)=nan;
end

best_scores=[scores(1,best_matches(1)), scores(2,best_matches(2)), scores(3,best_matches(3))];

end

