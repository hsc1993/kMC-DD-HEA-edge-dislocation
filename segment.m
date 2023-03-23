classdef segment

   properties
      idx_seg
      type
      rn_start
      rn_end
      idx_pairing_jog
   end

   methods
      function obj = MyClass(val)
         if nargin > 0
            obj.Prop = val;
         end
      end
   end
end