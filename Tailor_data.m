function ego_Tailored = Tailor_data(ego)
ego_Tailored=ego;

ego_ = fieldnames(ego);
for i=1:length(ego_)
    k = ego_(i);
    key = k{1};
    if length(ego.(key))==length(ego.t)
        Key=ego.(key);
        ego_Tailored.(key)=Key(1:end-50,1);
    end
   
end


end
