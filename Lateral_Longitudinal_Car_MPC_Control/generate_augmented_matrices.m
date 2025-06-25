function [A_aug, B_aug, C_aug, D_aug] = generate_augmented_matrices(Ad, Bd, Cd, Dd)

   A_aug = [Ad, Bd; zeros(size(Bd, 2), size(Ad, 2)), eye(size(Bd, 2))];
   B_aug = [Bd; eye(size(Bd, 2))];
   C_aug = [Cd, zeros(size(Cd, 1), size(Bd, 2))];
   D_aug = Dd;  % It is not used in the calculation
end