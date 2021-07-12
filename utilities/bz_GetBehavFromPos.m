function [behavior,trials] = bz_GetBehavFromPos(behavior, pos)

if size(pos,2) > 6
        behavior.position.x = pos(:,8);
        behavior.position.y = pos(:,10);
        behavior.position.z = pos(:,9);
        behavior.timestamps = pos(:,1);
        behavior.samplingRate = 1000 ./ mean(diff(behavior.timestamps))./1000;
        if nanstd(pos(:,7)) > 10  % determine unit of measure
            behavior.units = 'meters';
        else
            behavior.units = 'mm';
        end
        behavior.orientation.x = pos(:,4);
        behavior.orientation.y = pos(:,5);
        behavior.orientation.z = pos(:,6);
        behavior.orientation.w = pos(:,7);
        behavior.rotationType = 'quaternion';
        behavior.errorPerMarker = pos(:,11);
        
        c = ones(length(unique(behavior.events.trialConditions)),1);
        for i=1:size(behavior.events.trialIntervals,1)
            [a start] = min(abs(pos(:,1)-behavior.events.trialIntervals(i,1)));
            [a stop] = min(abs(pos(:,1)-behavior.events.trialIntervals(i,2)));
            
                behavior.events.trials{i}.x = pos(start:stop,8);
                behavior.events.trials{i}.y = pos(start:stop,10);
                behavior.events.trials{i}.z = pos(start:stop,9);
                behavior.events.trials{i}.orientation.x = pos(start:stop,4);
                behavior.events.trials{i}.orientation.y = pos(start:stop,5);
                behavior.events.trials{i}.orientation.z = pos(start:stop,6);
                behavior.events.trials{i}.orientation.w = pos(start:stop,7);
                behavior.events.trials{i}.errorPerMarker = pos(start:stop,11);
                behavior.rotationType = 'quaternion';
                behavior.events.trials{i}.timestamps = pos(start:stop,1);
                
            if stop-start > 10000
               error 
            end
            trials{behavior.events.trialConditions(i)}{c(behavior.events.trialConditions(i))} = pos(start:stop,:);
  
            c(behavior.events.trialConditions(i)) = 1+c(behavior.events.trialConditions(i));
        end
        
        [map mapping] = bz_trials2maps(trials);
        c = ones(length(unique(behavior.events.trialConditions)),1);
         for i=1:size(behavior.events.trialIntervals,1)
             behavior.events.map{behavior.events.trialConditions(i)}.x = map{behavior.events.trialConditions(i)}(:,8);
             behavior.events.map{behavior.events.trialConditions(i)}.y = map{behavior.events.trialConditions(i)}(:,10);
             behavior.events.map{behavior.events.trialConditions(i)}.z = map{behavior.events.trialConditions(i)}(:,9);
             behavior.events.trials{i}.mapping = mapping{behavior.events.trialConditions(i)}{c(behavior.events.trialConditions(i))}(:,13);
             
             c(behavior.events.trialConditions(i)) = 1+c(behavior.events.trialConditions(i));
         end
        
else
        behavior.position.x = mean(pos(:,[2 4]),2);
        behavior.position.y = mean(pos(:,[3 5]),2);
%         behavior.position.z = pos(:,9);
        behavior.timestamps = pos(:,1);
        behavior.samplingRate = 1000 ./ mean(diff(behavior.timestamps))./1000;
%         if nanstd(pos(:,7)) > 10  % determine unit of measure
%             behavior.units = 'meters';
%         else
%             behavior.units = 'mm';
%         end
%         behavior.orientation.x = pos(:,4);
%         behavior.orientation.y = pos(:,5);
%         behavior.orientation.z = pos(:,6);
%         behavior.orientation.w = pos(:,7);
%         behavior.rotationType = 'quaternion';
%         behavior.errorPerMarker = pos(:,11);
        
        c = ones(length(unique(behavior.events.trialConditions)),1);
        for i=1:size(behavior.events.trialIntervals,1)
            [a start] = min(abs(pos(:,1)-behavior.events.trialIntervals(i,1)));
            [a stop] = min(abs(pos(:,1)-behavior.events.trialIntervals(i,2)));
            
                behavior.events.trials{i}.x = mean(pos(start:stop,[2 4]),2);
                behavior.events.trials{i}.y = mean(pos(start:stop,[3 5]),2);
%                 behavior.events.trials{i}.z = pos(start:stop,9);
%                 behavior.events.trials{i}.orientation.x = pos(start:stop,4);
%                 behavior.events.trials{i}.orientation.y = pos(start:stop,5);
%                 behavior.events.trials{i}.orientation.z = pos(start:stop,6);
%                 behavior.events.trials{i}.orientation.w = pos(start:stop,7);
%                 behavior.events.trials{i}.errorPerMarker = pos(start:stop,11);
%                 behavior.rotationType = 'quaternion';
                behavior.events.trials{i}.timestamps = pos(start:stop,1);
                
            trials{behavior.events.trialConditions(i)}{c(behavior.events.trialConditions(i))} = pos(start:stop,:);
            
            c(behavior.events.trialConditions(i)) = 1+c(behavior.events.trialConditions(i));
        end
        [map mapping] = bz_trials2maps(trials);
        c = ones(length(unique(behavior.events.trialConditions)),1);
        for i=1:size(behavior.events.trialIntervals,1)
            behavior.events.map{behavior.events.trialConditions(i)}.x = map{behavior.events.trialConditions(i)}(:,8);
            behavior.events.map{behavior.events.trialConditions(i)}.y = map{behavior.events.trialConditions(i)}(:,10);
            behavior.events.map{behavior.events.trialConditions(i)}.z = map{behavior.events.trialConditions(i)}(:,9);
            behavior.events.trials{i}.mapping = mapping{behavior.events.trialConditions(i)}{c(behavior.events.trialConditions(i))}(:,13);

            c(behavior.events.trialConditions(i)) = 1+c(behavior.events.trialConditions(i));
        end
end


end