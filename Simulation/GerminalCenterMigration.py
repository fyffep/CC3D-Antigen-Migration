
from cc3d import CompuCellSetup
        


from GerminalCenterMigrationSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from GerminalCenterMigrationSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from GerminalCenterMigrationSteppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=100))




        
from GerminalCenterMigrationSteppables import UpdatePlotsSteppable
CompuCellSetup.register_steppable(steppable=UpdatePlotsSteppable(frequency=200))

CompuCellSetup.run()
