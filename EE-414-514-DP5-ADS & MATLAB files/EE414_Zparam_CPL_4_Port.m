function Z = EE414_ZParam_CPL_4_Port(Z0e, Z0o, theta_e, theta_o)

% Calculate the Z matrix elements
    Z11 = Z0e * cosh(theta_e) * cosh(theta_o) + Z0o * sinh(theta_e) * sinh(theta_o);
    Z12 = -Z0e * sinh(theta_e) * cosh(theta_o) + Z0o * cosh(theta_e) * sinh(theta_o);
    Z21 = -Z0e * sinh(theta_o) * cosh(theta_e) + Z0o * cosh(theta_o) * sinh(theta_e);
    Z22 = Z0e * cosh(theta_e) * cosh(theta_o) + Z0o * sinh(theta_e) * sinh(theta_o);

    % Construct the Z matrix
    Z = [Z11, Z12; Z21, Z22];
end
