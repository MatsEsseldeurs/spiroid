use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum WindModel {
    Disabled,
    #[default]
    Enabled,
}

impl WindModel {
    // Compute the magnetic torque if magnetism is enabled and the planet is inside the alfven radius.
    pub(crate) fn wind_torque(&self) -> bool {
        *self == Self::Enabled
    }
}
