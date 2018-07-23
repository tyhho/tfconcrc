function dspecies = deRepressionOED(t,species, parameters, I)

% Extract parameters from parameter array input
Vmax_P1 = parameters(1);
alpha_P2 = parameters(2);
Vmax_P2 = parameters(3);
Kr = parameters(4);
Ki = parameters(5);
n = parameters(6);
d_1 = parameters(7);
d_2 = parameters(8);

G = Ki^n / (Ki^n + I^n); %hill function for R-I interaction

dspecies = zeros(2,1);

dspecies(1) = Vmax_P1 - d_1 * species(1); %total R at steady state

dspecies(2) = alpha_P2 + Vmax_P2 * ( Kr /(Kr + species(1)*G) ) - d_2 * species(2);    %compute GFP output and store value

end