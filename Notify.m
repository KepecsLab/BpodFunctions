global BpodSystem

rig = BpodSystem.Data.Custom.Rig;
animal = BpodSystem.Data.Custom.Subject;
MailAlert({'Torben'},strcat(rig,':',animal,':PausedProtocol.'),'')