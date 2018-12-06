function group_isomorphism = lie_group_isomorphism (x)
% global center_position
global globalpara;
% center_radius = globalpara.center_radius;
radius = norm(x);
position = x;

cos_alpha = position(2) / radius;
sin_alpha = position(1) / radius;
group_isomorphism = [cos_alpha  -sin_alpha;
                     sin_alpha   cos_alpha];
end