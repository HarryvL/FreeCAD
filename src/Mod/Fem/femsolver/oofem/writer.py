# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2015 - Przemo Firszt <przemo@firszt.eu>                 *
# *   Copyright (c) 2015 - Bernd Hahnebach <bernd@bimstatik.org>            *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

__title__ = "OOFEM Writer"
__author__ = "Przemo Firszt, Bernd Hahnebach"
__url__ = "http://www.freecadweb.org"

## \addtogroup FEM
#  @{

import FreeCAD
import os
import sys
import time
import femmesh.meshtools as FemMeshTools
import feminout.importOOFEMMesh as importOOFEMMesh
from .. import writerbase as FemInputWriter
import six


class FemInputWriterOOFEM(FemInputWriter.FemInputWriter):
    def __init__(self,
                 analysis_obj, solver_obj,
                 mesh_obj, matlin_obj, matnonlin_obj,
                 fixed_obj, displacement_obj,
                 contact_obj, planerotation_obj, transform_obj,
                 selfweight_obj, force_obj, pressure_obj,
                 temperature_obj, heatflux_obj, initialtemperature_obj,
                 beamsection_obj, beamrotation_obj, shellthickness_obj, fluidsection_obj,
                 dir_name=None
                 ):

        FemInputWriter.FemInputWriter.__init__(
            self,
            analysis_obj, solver_obj,
            mesh_obj, matlin_obj, matnonlin_obj,
            fixed_obj, displacement_obj,
            contact_obj, planerotation_obj, transform_obj,
            selfweight_obj, force_obj, pressure_obj,
            temperature_obj, heatflux_obj, initialtemperature_obj,
            beamsection_obj, beamrotation_obj, shellthickness_obj, fluidsection_obj,
            dir_name)
        # self.dir_name does have a slash at the end
        self.main_file_name = self.mesh_object.Name + '.in'
        self.file_name = self.dir_name + self.main_file_name
        #self.FluidInletoutlet_ele = []
        #self.fluid_inout_nodes_file = self.dir_name + self.mesh_object.Name + '_inout_nodes.txt'
        FreeCAD.Console.PrintMessage('FemInputWriterCcx --> self.dir_name  -->  ' + self.dir_name + '\n')
        FreeCAD.Console.PrintMessage('FemInputWriterCcx --> self.main_file_name  -->  ' + self.main_file_name + '\n')
        FreeCAD.Console.PrintMessage('FemInputWriterCcx --> self.file_name  -->  ' + self.file_name + '\n')

    def write_OOFEM_input_file(self):
        
        timestart = time.clock()
        
        inpfile = open(self.file_name, 'a')        
        
        self.write_output_file_record(inpfile)
        ''' Filename String '''
        
        self.write_job_description_record(inpfile)
        ''' Job description string '''
        
        self.write_analysis_record(inpfile)
        '''*AnalysisType     nsteps # (in)
                            [renumber # (in) ]
                            [profileopt # (in) ]
                            attributes # (string)
                            [ninitmodules # (in) ]
                            [nmodules # (in) ]
                            [nxfemman # (in) ]
        
                            *StaticStructural
                            nsteps # (in)
                            [deltat # (...) ]
                            [prescribedtimes # (...) ]
                            [stiffmode # (...) ]
                            [nonlocalext # (...) ]
                            [sparselinsolverparams # (...) ]
                            
                            *LinearStability
                             nroot # (in)
                             rtolv # (rn)
                             [eigensolverparams # (...) ]
                             
                            *NonLinearStatic
                             [nmsteps # (in) ]
                             nsteps # (in)
                             [contextOutputStep # (in) ]
                             [sparselinsolverparams # (string) ]
                             [nonlinform # (in) ]
                             [nonlocstiff # (in) ]
                             [nonlocalext]
                             [loadbalancing]'''

        self.write_domain_record(inpfile)

        '''domain *domainType
                              *2dPlaneStress
                              *2d-Truss
                              *3d
                              *2dMindlinPlate
                              *3dShell
                              *2dBeam
        
        '''
        
        self.write_output_manager_record(inpfile)
        
        '''OutputManager
           [tstep all]
           [tstep step # (in) ]
           [tsteps out # (rl) ]
           [dofman all]
           [dofman output # (rl) ]
           [dofman except # (rl) ]
           [element all]
           [element output # (rl) ]
           [element except # (rl) ]
           ndofman # (in)
           nelem # (in)
           ncrosssect # (in)
           nmat # (in)
           nbc # (in)
           nic # (in)
           nltf # (in)
           [nbarrier # (in) ]'''
           
        self.write_dof_manager_record(inpfile)
           
        '''*DofManagerType
            (num#) (in)
            [load # (ra) ]
            [DofIDMask # (ia) ]
            [bc # (ia) ]
            [ic # (ia) ]
            [doftype # (ia) masterMask # (ia) ]
            [shared]i | h[remote]i | h[null]
            [partitions # (ia) ]
                                                 *Node coords # (ra)
                                                  [lcs # (ra) ]'''
                                                  
        self.write_element_record(inpfile)
           
        '''*ElementType
            (num#) (in)
            mat # (in) crossSect # (in) nodes # (ia)
            [bodyLoads # (ia) ] [boundaryLoads # (ia) ]
            [activityltf # (in) ] [lcs # (ra) ]
            [partitions # (ia) ] [remote]
                                                 *beam2d
                                                  2D beam element
                                                  [dofstocondense # (ia) ]
                                                 *beam3d
                                                  3D beam element
                                                  refnode # (in) [dofstocondense # (ia) ]
                                                 *planestress2d
                                                  4-noded 2D quadrilateral element for plane stress analysis
                                                  [NIP # (in) ]
                                                 *qplanestress2d
                                                  8-noded 2D quadrilateral element for plane stress analysis
                                                  [NIP # (in) ]
                                                 *trplanestress2d
                                                  3-noded 2D triangular element for plane stress analysis                                                 
                                                 *qtrplstr
                                                  6-noded 2D triangular element for plane stress analysis [NIP # (in) ]
                                                 *quad1planestrain
                                                  4-noded 2D quadrilateral element for plane strain analysis
                                                  [NIP # (in) ]
                                                 *trplanestrain
                                                  3-noded 2D triangular element for plane strain analysis
                                                 *lspace
                                                  Linear 8-node isoparametric brick element
                                                  [NIP # (in) ]
                                                 *qspace
                                                  Quadratic 20-node isoparametric brick element
                                                  [NIP # (in) ]
                                                 *LTRSpace
                                                  Linear 4-node tetrahedra element
                                                 *QTRSpace
                                                  Quadratic 10-node tetrahedra element
                                                  [NIP # (in) ] '''

        self.write_set_record(inpfile)
        
        '''Set (num#) (in)
           [elements # (ia) ] [elementranges # (rl) ] [allElements]
           [nodes # (ia) ] [noderanges # (rl) ] [allNodes]
           [elementboundaries # (ia) ] [elementedges # (ia) ] '''
           
        self.write_cross_section_record(inpfile)
        
        '''*CrossSectType (num#) (in)
                                        *SimpleCS [thick # (rn) ] [width # (rn) ] [area # (rn) ]
                                         [iy # (rn) ] [iz # (rn) ] [ik # (rn) ]
                                         [shearareay # (rn) ] [shearareaz # (rn) ] beamshearcoeff # (rn)
                                         
                                        *VariableCS [thick # (expr) ] [width # (expr) ] [area # (expr) ]
                                        [iy # (expr) ] [iz # (expr) ] [ik # (expr) ]
                                        [shearareay # (expr) ] [shearareaz # (expr) ]
                                        
                                        *LayeredCS nLayers # (in)
                                        LayerMaterials # (ia)
                                        Thicks # (ra) Widths # (ra)
                                        midSurf # (rn) '''

        self.write_material_type_record(inpfile)
        
        '''*MaterialType (num#) (in) d # (rn)
                                         Linear isotropic elastic material
                                        *IsoLE num (in) # d (rn) # E (rn) # n (rn) # tAlpha (rn) #
                                         Mooney-Rivlin
                                        *MooneyRivlin (in) # d (rn) # K (rn) # C1 (rn) # C2 (rn) #
                                         Large-strain master material material
                                        *LSmasterMat (in) # m (rn) # slavemat (in) #
                                         DP material
                                        *DruckerPrager num (in) # d (rn) # tAlpha (rn) # E (rn) # n (rn) # alpha (rn) # alphaPsi (rn) # ht (in) # iys (rn) # lys (rn) # hm (rn) # kc (rn) # [ yieldtol (rn) #]
                                         Mises plasticity model with isotropic hardening
                                        *MisesMat (in) # d (rn) # E (rn) # n (rn) # sig0 (rn) # H (rn) # omega crit (rn) #a (rn) #
                                         Rotating crack model for concrete
                                        *Concrete3 d (rn) # E (rn) # n (rn) # Gf (rn) # Ft (rn) # exp soft (in) # tAlpha (rn) #
                                         EC2CreepMat model for concrete creep and shrinkage
                                        *EC2CreepMat n (rn) # [ begOfTimeOfInterest (rn) #] [ end-OfTimeOfInterest (rn) #] relMatAge (rn) # [ timeFactor (rn) #] stiffnessFactor (rn) # [ tAlpha (rn) #] fcm28 (rn) # t0 (rn) # cem- Type (in) # [ henv (rn) #] h0 (rn) # shType (in) # [ spectrum ][ temperatureDependent ]'''
           
                   
        # footer
        self.write_footer(inpfile)
        inpfile.close()
        FreeCAD.Console.PrintMessage("Writing time input file: " + str(time.clock() - timestart) + ' \n\n')
        return self.file_name


    def write_output_file_record(self, f):
        f.write('Filename String\n')

    def write_job_description_record(self, f):
        f.write('Job description string\n')

    def write_analysis_record(self, f):
        f.write('#\n')
        f.write('# Analysis Record\n')
        f.write('#\n')

    def write_domain_record(self, f):
        f.write('#\n')
        f.write('# Domain Record\n')
        f.write('#\n')

    def write_output_manager_record(self, f):
        f.write('#\n')
        f.write('# Output Manager Record\n')
        f.write('#\n')

    def write_dof_manager_record(self, f):
        f.write('#\n')
        f.write('# DOF Manager Record\n')
        f.write('#\n')

    def write_element_record(self, f):
        self.OOFEM_element_type=16
        f.write('#\n')
        f.write('# Element Record\n')
        f.write('#\n')
        f.write("self.femnodes_mesh {0}\n".format(self.femnodes_mesh))
        f.write("self.femelement_table {0}\n".format(self.femelement_table))
        importOOFEMMesh.write_OOFEM_mesh_to_file(self.femnodes_mesh, self.femelement_table, self.OOFEM_element_type, f)
        
    def write_set_record(self, f):
        f.write('#\n')
        f.write('# Set Record\n')
        f.write('#\n')

    def write_cross_section_record(self, f):
        f.write('#\n')
        f.write('# Cross Section Record\n')
        f.write('#\n')

    def write_material_type_record(self, f):
        f.write('#\n')
        f.write('# Material Type Record\n')
        f.write('#\n')


    def write_footer(self, f):
        f.write('#\n')
        f.write('# End of Input File\n')
        f.write('#\n')



    def write_element_sets_material_and_femelement_type(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Element sets for materials and FEM element type (solid, shell, beam, fluid)\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))

        # in any case if we have beams, we're going to need the element ids for the rotation elsets
        if self.beamsection_objects:
            # we will need to split the beam even for one beamobj
            # because no beam in z-direction can be used in ccx without a special adjustment
            # thus they need an own ccx_elset
            self.get_element_rotation1D_elements()

        # get the element ids for face and edge elements and write them into the objects
        if len(self.shellthickness_objects) > 1:
            self.get_element_geometry2D_elements()
        if len(self.beamsection_objects) > 1:
            self.get_element_geometry1D_elements()
        if len(self.fluidsection_objects) > 1:
            self.get_element_fluid1D_elements()

        # get the element ids for material objects and write them into the material object
        if len(self.material_objects) > 1:
            self.get_material_elements()

        # create the ccx_elsets
        if len(self.material_objects) == 1:
            if self.femmesh.Volumes:
                # we only could do this for volumes, if a mesh contains volumes
                # we're going to use them in the analysis
                # but a mesh could contain the element faces of the volumes as faces
                # and the edges of the faces as edges, there we have to check for some geometric objects
                self.get_ccx_elsets_single_mat_solid()
            if len(self.shellthickness_objects) == 1:
                self.get_ccx_elsets_single_mat_single_shell()
            elif len(self.shellthickness_objects) > 1:
                self.get_ccx_elsets_single_mat_multiple_shell()
            if len(self.beamsection_objects) == 1:
                self.get_ccx_elsets_single_mat_single_beam()
            elif len(self.beamsection_objects) > 1:
                self.get_ccx_elsets_single_mat_multiple_beam()
            if len(self.fluidsection_objects) == 1:
                self.get_ccx_elsets_single_mat_single_fluid()
            elif len(self.fluidsection_objects) > 1:
                self.get_ccx_elsets_single_mat_multiple_fluid()
        elif len(self.material_objects) > 1:
            if self.femmesh.Volumes:
                # we only could do this for volumes, if a mseh contains volumes
                # we're going to use them in the analysis
                # but a mesh could contain the element faces of the volumes as faces
                # and the edges of the faces as edges, there we have to check for some geometric objects
                self.get_ccx_elsets_multiple_mat_solid()  # volume is a bit special, because retrieving ids from group mesh data is implemented
            if len(self.shellthickness_objects) == 1:
                self.get_ccx_elsets_multiple_mat_single_shell()
            elif len(self.shellthickness_objects) > 1:
                self.get_ccx_elsets_multiple_mat_multiple_shell()
            if len(self.beamsection_objects) == 1:
                self.get_ccx_elsets_multiple_mat_single_beam()
            elif len(self.beamsection_objects) > 1:
                self.get_ccx_elsets_multiple_mat_multiple_beam()
            if len(self.fluidsection_objects) == 1:
                self.get_ccx_elsets_multiple_mat_single_fluid()
            elif len(self.fluidsection_objects) > 1:
                self.get_ccx_elsets_multiple_mat_multiple_fluid()

        # TODO: some elemetIDs are collected for 1D-Flow calculation,
        # this should be a def somewhere else, preferable inside the get_ccx_elsets_... methods
        for ccx_elset in self.ccx_elsets:
            if ccx_elset['ccx_elset'] and not isinstance(ccx_elset['ccx_elset'], six.string_types):  # use six to be sure to be Python 2.7 and 3.x compatible
                if 'fluidsection_obj'in ccx_elset:
                    fluidsec_obj = ccx_elset['fluidsection_obj']
                    if fluidsec_obj.SectionType == 'Liquid':
                        if (fluidsec_obj.LiquidSectionType == "PIPE INLET") or (fluidsec_obj.LiquidSectionType == "PIPE OUTLET"):
                            elsetchanged = False
                            counter = 0
                            for elid in ccx_elset['ccx_elset']:
                                counter = counter + 1
                                if (elsetchanged is False) and (fluidsec_obj.LiquidSectionType == "PIPE INLET"):
                                    self.FluidInletoutlet_ele.append([str(elid), fluidsec_obj.LiquidSectionType, 0])  # 3rd index is to track which line number the element is defined
                                    elsetchanged = True
                                elif (fluidsec_obj.LiquidSectionType == "PIPE OUTLET") and (counter == len(ccx_elset['ccx_elset'])):
                                    self.FluidInletoutlet_ele.append([str(elid), fluidsec_obj.LiquidSectionType, 0])  # 3rd index is to track which line number the element is defined

        # write ccx_elsets to file
        for ccx_elset in self.ccx_elsets:
            f.write('*ELSET,ELSET=' + ccx_elset['ccx_elset_name'] + '\n')
            if isinstance(ccx_elset['ccx_elset'], six.string_types):  # use six to be sure to be Python 2.7 and 3.x compatible
                f.write(ccx_elset['ccx_elset'] + '\n')
            else:
                for elid in ccx_elset['ccx_elset']:
                    f.write(str(elid) + ',\n')

    def write_node_sets_constraints_fixed(self, f):
        # get nodes
        self.get_constraints_fixed_nodes()
        # write nodes to file
        f.write('\n***********************************************************\n')
        f.write('** Node sets for fixed constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.fixed_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            fix_obj = femobj['Object']
            f.write('** ' + fix_obj.Label + '\n')
            if self.femmesh.Volumes and (len(self.shellthickness_objects) > 0 or len(self.beamsection_objects) > 0):
                if len(femobj['NodesSolid']) > 0:
                    f.write('*NSET,NSET=' + fix_obj.Name + 'Solid\n')
                    for n in femobj['NodesSolid']:
                        f.write(str(n) + ',\n')
                if len(femobj['NodesFaceEdge']) > 0:
                    f.write('*NSET,NSET=' + fix_obj.Name + 'FaceEdge\n')
                    for n in femobj['NodesFaceEdge']:
                        f.write(str(n) + ',\n')
            else:
                f.write('*NSET,NSET=' + fix_obj.Name + '\n')
                for n in femobj['Nodes']:
                    f.write(str(n) + ',\n')

    def write_node_sets_constraints_displacement(self, f):
        # get nodes
        self.get_constraints_displacement_nodes()
        # write nodes to file
        f.write('\n***********************************************************\n')
        f.write('** Node sets for prescribed displacement constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.displacement_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            disp_obj = femobj['Object']
            f.write('** ' + disp_obj.Label + '\n')
            f.write('*NSET,NSET=' + disp_obj.Name + '\n')
            for n in femobj['Nodes']:
                f.write(str(n) + ',\n')

    def write_node_sets_constraints_planerotation(self, f):
        # get nodes
        self.get_constraints_planerotation_nodes()
        # write nodes to file
        if not self.femnodes_mesh:
            self.femnodes_mesh = self.femmesh.Nodes
        f.write('\n***********************************************************\n')
        f.write('** Node sets for plane rotation constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        # info about self.constraint_conflict_nodes:
        # is used to check if MPC and constraint fixed and constraint displacement share same nodes,
        # because MPC's and constriants fixed and constraints displacement can't share same nodes.
        # Thus call write_node_sets_constraints_planerotation has to be after constraint fixed and constraint displacement
        for femobj in self.planerotation_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            l_nodes = femobj['Nodes']
            fric_obj = femobj['Object']
            f.write('** ' + fric_obj.Label + '\n')
            f.write('*NSET,NSET=' + fric_obj.Name + '\n')
            # Code to extract nodes and coordinates on the PlaneRotation support face
            nodes_coords = []
            for node in l_nodes:
                nodes_coords.append((node, self.femnodes_mesh[node].x, self.femnodes_mesh[node].y, self.femnodes_mesh[node].z))
            node_planerotation = FemMeshTools.get_three_non_colinear_nodes(nodes_coords)
            for i in range(len(l_nodes)):
                if l_nodes[i] not in node_planerotation:
                    node_planerotation.append(l_nodes[i])
            MPC_nodes = []
            for i in range(len(node_planerotation)):
                cnt = 0
                for j in range(len(self.constraint_conflict_nodes)):
                    if node_planerotation[i] == self.constraint_conflict_nodes[j]:
                        cnt = cnt + 1
                if cnt == 0:
                    MPC = node_planerotation[i]
                    MPC_nodes.append(MPC)
            for i in range(len(MPC_nodes)):
                f.write(str(MPC_nodes[i]) + ',\n')

    def write_surfaces_contraints_contact(self, f):
        # get surface nodes and write them to file
        f.write('\n***********************************************************\n')
        f.write('** Surfaces for contact constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        obj = 0
        for femobj in self.contact_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            contact_obj = femobj['Object']
            f.write('** ' + contact_obj.Label + '\n')
            cnt = 0
            obj = obj + 1
            for o, elem_tup in contact_obj.References:
                for elem in elem_tup:
                    ref_shape = o.Shape.getElement(elem)
                    cnt = cnt + 1
                    if ref_shape.ShapeType == 'Face':
                        if cnt == 1:
                            name = "DEP" + str(obj)
                        else:
                            name = "IND" + str(obj)
                        f.write('*SURFACE, NAME =' + name + '\n')
                        v = self.mesh_object.FemMesh.getccxVolumesByFace(ref_shape)
                        for i in v:
                            f.write("{},S{}\n".format(i[0], i[1]))

    def write_node_sets_constraints_transform(self, f):
        # get nodes
        self.get_constraints_transform_nodes()
        # write nodes to file
        f.write('\n***********************************************************\n')
        f.write('** Node sets for transform constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.transform_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            trans_obj = femobj['Object']
            f.write('** ' + trans_obj.Label + '\n')
            if trans_obj.TransformType == "Rectangular":
                f.write('*NSET,NSET=Rect' + trans_obj.Name + '\n')
            elif trans_obj.TransformType == "Cylindrical":
                f.write('*NSET,NSET=Cylin' + trans_obj.Name + '\n')
            for n in femobj['Nodes']:
                f.write(str(n) + ',\n')

    def write_node_sets_constraints_temperature(self, f):
        # get nodes
        self.get_constraints_temperature_nodes()
        # write nodes to file
        f.write('\n***********************************************************\n')
        f.write('** Node sets for temperature constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.temperature_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            temp_obj = femobj['Object']
            f.write('** ' + temp_obj.Label + '\n')
            f.write('*NSET,NSET=' + temp_obj.Name + '\n')
            for n in femobj['Nodes']:
                f.write(str(n) + ',\n')

    def write_materials(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Materials\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        f.write('** Young\'s modulus unit is MPa = N/mm2\n')
        if self.analysis_type == "frequency" or self.selfweight_objects or (self.analysis_type == "thermomech" and not self.solver_obj.ThermoMechSteadyState):
            f.write('** Density\'s unit is t/mm^3\n')
        if self.analysis_type == "thermomech":
            f.write('** Thermal conductivity unit is kW/mm/K = t*mm/K*s^3\n')
            f.write('** Specific Heat unit is kJ/t/K = mm^2/s^2/K\n')
        for femobj in self.material_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            mat_obj = femobj['Object']
            mat_info_name = mat_obj.Material['Name']
            mat_name = mat_obj.Name
            mat_label = mat_obj.Label
            # get material properties of solid material, Currently in SI units: M/kg/s/Kelvin
            if mat_obj.Category == 'Solid':
                YM = FreeCAD.Units.Quantity(mat_obj.Material['YoungsModulus'])
                YM_in_MPa = float(YM.getValueAs('MPa'))
                PR = float(mat_obj.Material['PoissonRatio'])
            if self.analysis_type == "frequency" or self.selfweight_objects or (self.analysis_type == "thermomech" and not self.solver_obj.ThermoMechSteadyState):
                density = FreeCAD.Units.Quantity(mat_obj.Material['Density'])
                density_in_tonne_per_mm3 = float(density.getValueAs('t/mm^3'))
            if self.analysis_type == "thermomech":
                TC = FreeCAD.Units.Quantity(mat_obj.Material['ThermalConductivity'])
                TC_in_WmK = float(TC.getValueAs('W/m/K'))  # SvdW: Add factor to force units to results' base units of t/mm/s/K - W/m/K results in no factor needed
                SH = FreeCAD.Units.Quantity(mat_obj.Material['SpecificHeat'])
                SH_in_JkgK = float(SH.getValueAs('J/kg/K')) * 1e+06  # SvdW: Add factor to force units to results' base units of t/mm/s/K
                if mat_obj.Category == 'Solid':
                    TEC = FreeCAD.Units.Quantity(mat_obj.Material['ThermalExpansionCoefficient'])
                    TEC_in_mmK = float(TEC.getValueAs('mm/mm/K'))
                elif mat_obj.Category == 'Fluid':
                    DV = FreeCAD.Units.Quantity(mat_obj.Material['DynamicViscosity'])
                    DV_in_tmms = float(DV.getValueAs('t/mm/s'))
            # write material properties
            f.write('** FreeCAD material name: ' + mat_info_name + '\n')
            f.write('** ' + mat_label + '\n')
            f.write('*MATERIAL, NAME=' + mat_name + '\n')
            if mat_obj.Category == 'Solid':
                f.write('*ELASTIC\n')
                f.write('{0:.0f}, {1:.3f}\n'.format(YM_in_MPa, PR))

            if self.analysis_type == "frequency" or self.selfweight_objects or (self.analysis_type == "thermomech" and not self.solver_obj.ThermoMechSteadyState):
                f.write('*DENSITY\n')
                f.write('{0:.3e}\n'.format(density_in_tonne_per_mm3))
            if self.analysis_type == "thermomech":
                if mat_obj.Category == 'Solid':
                    f.write('*CONDUCTIVITY\n')
                    f.write('{0:.3f}\n'.format(TC_in_WmK))
                    f.write('*EXPANSION\n')
                    f.write('{0:.3e}\n'.format(TEC_in_mmK))
                    f.write('*SPECIFIC HEAT\n')
                    f.write('{0:.3e}\n'.format(SH_in_JkgK))
                elif mat_obj.Category == 'Fluid':
                    f.write('*FLUID CONSTANTS\n')
                    f.write('{0:.3e}, {1:.3e}\n'.format(SH_in_JkgK, DV_in_tmms))

            # nonlinear material properties
            if self.solver_obj.MaterialNonlinearity == 'nonlinear':
                for femobj in self.material_nonlinear_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
                    nl_mat_obj = femobj['Object']
                    if nl_mat_obj.LinearBaseMaterial == mat_obj:
                        if nl_mat_obj.MaterialModelNonlinearity == "simple hardening":
                            f.write('*PLASTIC\n')
                            if nl_mat_obj.YieldPoint1:
                                f.write(nl_mat_obj.YieldPoint1 + '\n')
                            if nl_mat_obj.YieldPoint2:
                                f.write(nl_mat_obj.YieldPoint2 + '\n')
                            if nl_mat_obj.YieldPoint3:
                                f.write(nl_mat_obj.YieldPoint3 + '\n')
                    f.write('\n')

    def write_constraints_initialtemperature(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Initial temperature constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        f.write('*INITIAL CONDITIONS,TYPE=TEMPERATURE\n')
        for itobj in self.initialtemperature_objects:  # Should only be one
            inittemp_obj = itobj['Object']
            f.write('{0},{1}\n'.format(self.ccx_nall, inittemp_obj.initialTemperature))  # OvG: Initial temperature

    def write_femelementsets(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Sections\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for ccx_elset in self.ccx_elsets:
            if ccx_elset['ccx_elset']:
                if 'beamsection_obj'in ccx_elset:  # beam mesh
                    beamsec_obj = ccx_elset['beamsection_obj']
                    elsetdef = 'ELSET=' + ccx_elset['ccx_elset_name'] + ', '
                    material = 'MATERIAL=' + ccx_elset['mat_obj_name']
                    normal = ccx_elset['beam_normal']
                    if beamsec_obj.SectionType == 'Rectangular':
                        height = beamsec_obj.RectHeight.getValueAs('mm')
                        width = beamsec_obj.RectWidth.getValueAs('mm')
                        section_type = ', SECTION=RECT'
                        setion_geo = str(height) + ', ' + str(width) + '\n'
                        setion_def = '*BEAM SECTION, ' + elsetdef + material + section_type + '\n'
                        setion_nor = str(normal[0]) + ', ' + str(normal[1]) + ', ' + str(normal[2]) + '\n'
                    elif beamsec_obj.SectionType == 'Circular':
                        radius = 0.5 * beamsec_obj.CircDiameter.getValueAs('mm')
                        section_type = ', SECTION=CIRC'
                        setion_geo = str(radius) + '\n'
                        setion_def = '*BEAM SECTION, ' + elsetdef + material + section_type + '\n'
                        setion_nor = str(normal[0]) + ', ' + str(normal[1]) + ', ' + str(normal[2]) + '\n'
                    elif beamsec_obj.SectionType == 'Pipe':
                        radius = 0.5 * beamsec_obj.PipeDiameter.getValueAs('mm')
                        thickness = beamsec_obj.PipeThickness.getValueAs('mm')
                        section_type = ', SECTION=PIPE'
                        setion_geo = str(radius) + ', ' + str(thickness) + '\n'
                        setion_def = '*BEAM GENERAL SECTION, ' + elsetdef + material + section_type + '\n'
                        setion_nor = str(normal[0]) + ', ' + str(normal[1]) + ', ' + str(normal[2]) + '\n'
                    f.write(setion_def)
                    f.write(setion_geo)
                    f.write(setion_nor)
                elif 'fluidsection_obj'in ccx_elset:  # fluid mesh
                    fluidsec_obj = ccx_elset['fluidsection_obj']
                    elsetdef = 'ELSET=' + ccx_elset['ccx_elset_name'] + ', '
                    material = 'MATERIAL=' + ccx_elset['mat_obj_name']
                    if fluidsec_obj.SectionType == 'Liquid':
                        section_type = fluidsec_obj.LiquidSectionType
                        if (section_type == "PIPE INLET") or (section_type == "PIPE OUTLET"):
                            section_type = "PIPE INOUT"
                        setion_def = '*FLUID SECTION, ' + elsetdef + 'TYPE=' + section_type + ', ' + material + '\n'
                        setion_geo = liquid_section_def(fluidsec_obj, section_type)
                    elif fluidsec_obj.SectionType == 'Gas':
                        section_type = fluidsec_obj.GasSectionType
                    elif fluidsec_obj.SectionType == 'Open Channel':
                        section_type = fluidsec_obj.ChannelSectionType
                    f.write(setion_def)
                    f.write(setion_geo)
                elif 'shellthickness_obj'in ccx_elset:  # shell mesh
                    shellth_obj = ccx_elset['shellthickness_obj']
                    elsetdef = 'ELSET=' + ccx_elset['ccx_elset_name'] + ', '
                    material = 'MATERIAL=' + ccx_elset['mat_obj_name']
                    setion_def = '*SHELL SECTION, ' + elsetdef + material + '\n'
                    setion_geo = str(shellth_obj.Thickness.getValueAs('mm')) + '\n'
                    f.write(setion_def)
                    f.write(setion_geo)
                else:  # solid mesh
                    elsetdef = 'ELSET=' + ccx_elset['ccx_elset_name'] + ', '
                    material = 'MATERIAL=' + ccx_elset['mat_obj_name']
                    setion_def = '*SOLID SECTION, ' + elsetdef + material + '\n'
                    f.write(setion_def)

    def write_step_begin(self, f):
        f.write('\n***********************************************************\n')
        f.write('** At least one step is needed to run an CalculiX analysis of FreeCAD\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        # STEP line
        step = '*STEP'
        if self.solver_obj.GeometricalNonlinearity == "nonlinear":
            if self.analysis_type == 'static' or self.analysis_type == 'thermomech':
                step += ', NLGEOM'   # https://www.comsol.com/blogs/what-is-geometric-nonlinearity/
            elif self.analysis_type == 'frequency':
                FreeCAD.Console.PrintMessage('Analysis type frequency and geometrical nonlinear analysis are not allowed together, linear is used instead!\n')
        if self.solver_obj.IterationsThermoMechMaximum:
            if self.analysis_type == 'thermomech':
                step += ', INC=' + str(self.solver_obj.IterationsThermoMechMaximum)
            elif self.analysis_type == 'static' or self.analysis_type == 'frequency':
                pass  # parameter is for thermomechanical analysis only, see ccx manual *STEP
        # write step line
        f.write(step + '\n')
        # CONTROLS line
        # all analysis types, ... really in frequency too?!?
        if self.solver_obj.IterationsControlParameterTimeUse:
            f.write('*CONTROLS, PARAMETERS=TIME INCREMENTATION\n')
            f.write(self.solver_obj.IterationsControlParameterIter + '\n')
            f.write(self.solver_obj.IterationsControlParameterCutb + '\n')
        # ANALYSIS type line
        # analysis line --> analysis type
        if self.analysis_type == 'static':
            analysis_type = '*STATIC'
        elif self.analysis_type == 'frequency':
            analysis_type = '*FREQUENCY'
        elif self.analysis_type == 'thermomech':
            analysis_type = '*COUPLED TEMPERATURE-DISPLACEMENT'
        # analysis line --> solver type
        if self.solver_obj.MatrixSolverType == "default":
            pass
        elif self.solver_obj.MatrixSolverType == "spooles":
            analysis_type += ', SOLVER=SPOOLES'
        elif self.solver_obj.MatrixSolverType == "iterativescaling":
            analysis_type += ', SOLVER=ITERATIVE SCALING'
        elif self.solver_obj.MatrixSolverType == "iterativecholesky":
            analysis_type += ', SOLVER=ITERATIVE CHOLESKY'
        # analysis line --> user defined incrementations --> parameter DIRECT --> completely switch off ccx automatic incrementation
        if self.solver_obj.IterationsUserDefinedIncrementations:
            if self.analysis_type == 'static':
                analysis_type += ', DIRECT'
            elif self.analysis_type == 'thermomech':
                analysis_type += ', DIRECT'
            elif self.analysis_type == 'frequency':
                FreeCAD.Console.PrintMessage('Analysis type frequency and IterationsUserDefinedIncrementations are not allowed together, it is ignored\n')
        # analysis line --> steadystate --> thermomech only
        if self.solver_obj.ThermoMechSteadyState:
            if self.analysis_type == 'thermomech':  # bernd: I do not know if STEADY STATE is allowed with DIRECT but since time steps are 1.0 it makes no sense IMHO
                analysis_type += ', STEADY STATE'
                self.solver_obj.TimeInitialStep = 1.0  # Set time to 1 and ignore user inputs for steady state
                self.solver_obj.TimeEnd = 1.0
            elif self.analysis_type == 'static' or self.analysis_type == 'frequency':
                pass  # not supported for static and frequency!
        # ANALYSIS parameter line
        analysis_parameter = ''
        if self.analysis_type == 'static':
            if self.solver_obj.IterationsUserDefinedIncrementations is True or self.solver_obj.IterationsUserDefinedTimeStepLength is True:
                analysis_parameter = '{},{}'.format(self.solver_obj.TimeInitialStep, self.solver_obj.TimeEnd)
        elif self.analysis_type == 'frequency':
            if self.solver_obj.EigenmodeLowLimit == 0.0 and self.solver_obj.EigenmodeHighLimit == 0.0:
                analysis_parameter = '{}\n'.format(self.solver_obj.EigenmodesCount)
            else:
                analysis_parameter = '{},{},{}\n'.format(self.solver_obj.EigenmodesCount, self.solver_obj.EigenmodeLowLimit, self.solver_obj.EigenmodeHighLimit)
        elif self.analysis_type == 'thermomech':
            analysis_parameter = '{},{}'.format(self.solver_obj.TimeInitialStep, self.solver_obj.TimeEnd)  # OvG: 1.0 increment, total time 1 for steady state will cut back automatically
        # write analysis type line, analysis parameter line
        f.write(analysis_type + '\n')
        f.write(analysis_parameter + '\n')

    def write_constraints_fixed(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Fixed Constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.fixed_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            f.write('** ' + femobj['Object'].Label + '\n')
            fix_obj_name = femobj['Object'].Name
            if self.femmesh.Volumes and (len(self.shellthickness_objects) > 0 or len(self.beamsection_objects) > 0):
                if len(femobj['NodesSolid']) > 0:
                    f.write('*BOUNDARY\n')
                    f.write(fix_obj_name + 'Solid' + ',1\n')
                    f.write(fix_obj_name + 'Solid' + ',2\n')
                    f.write(fix_obj_name + 'Solid' + ',3\n')
                    f.write('\n')
                if len(femobj['NodesFaceEdge']) > 0:
                    f.write('*BOUNDARY\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',1\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',2\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',3\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',4\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',5\n')
                    f.write(fix_obj_name + 'FaceEdge' + ',6\n')
                    f.write('\n')
            else:
                f.write('*BOUNDARY\n')
                f.write(fix_obj_name + ',1\n')
                f.write(fix_obj_name + ',2\n')
                f.write(fix_obj_name + ',3\n')
                if self.beamsection_objects or self.shellthickness_objects:
                    f.write(fix_obj_name + ',4\n')
                    f.write(fix_obj_name + ',5\n')
                    f.write(fix_obj_name + ',6\n')
                f.write('\n')

    def write_constraints_displacement(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Displacement constraint applied\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.displacement_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            f.write('** ' + femobj['Object'].Label + '\n')
            disp_obj = femobj['Object']
            disp_obj_name = disp_obj.Name
            f.write('*BOUNDARY\n')
            if disp_obj.xFix:
                f.write(disp_obj_name + ',1\n')
            elif not disp_obj.xFree:
                f.write(disp_obj_name + ',1,1,' + str(disp_obj.xDisplacement) + '\n')
            if disp_obj.yFix:
                f.write(disp_obj_name + ',2\n')
            elif not disp_obj.yFree:
                f.write(disp_obj_name + ',2,2,' + str(disp_obj.yDisplacement) + '\n')
            if disp_obj.zFix:
                f.write(disp_obj_name + ',3\n')
            elif not disp_obj.zFree:
                f.write(disp_obj_name + ',3,3,' + str(disp_obj.zDisplacement) + '\n')

            if self.beamsection_objects or self.shellthickness_objects:
                if disp_obj.rotxFix:
                    f.write(disp_obj_name + ',4\n')
                elif not disp_obj.rotxFree:
                    f.write(disp_obj_name + ',4,4,' + str(disp_obj.xRotation) + '\n')
                if disp_obj.rotyFix:
                    f.write(disp_obj_name + ',5\n')
                elif not disp_obj.rotyFree:
                    f.write(disp_obj_name + ',5,5,' + str(disp_obj.yRotation) + '\n')
                if disp_obj.rotzFix:
                    f.write(disp_obj_name + ',6\n')
                elif not disp_obj.rotzFree:
                    f.write(disp_obj_name + ',6,6,' + str(disp_obj.zRotation) + '\n')
        f.write('\n')

    def write_constraints_contact(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Contact Constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        obj = 0
        for femobj in self.contact_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            obj = obj + 1
            contact_obj = femobj['Object']
            f.write('** ' + contact_obj.Label + '\n')
            f.write('*CONTACT PAIR, INTERACTION=INT' + str(obj) + ',TYPE=SURFACE TO SURFACE\n')
            ind_surf = "IND" + str(obj)
            dep_surf = "DEP" + str(obj)
            f.write(dep_surf + ',' + ind_surf + '\n')
            f.write('*SURFACE INTERACTION, NAME=INT' + str(obj) + '\n')
            f.write('*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=LINEAR\n')
            slope = contact_obj.Slope
            f.write(str(slope) + ' \n')
            friction = contact_obj.Friction
            if friction > 0:
                f.write('*FRICTION \n')
                stick = (slope / 10.0)
                f.write(str(friction) + ', ' + str(stick) + ' \n')

    def write_constraints_planerotation(self, f):
        f.write('\n***********************************************************\n')
        f.write('** PlaneRotation Constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.planerotation_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            f.write('** ' + femobj['Object'].Label + '\n')
            fric_obj_name = femobj['Object'].Name
            f.write('*MPC\n')
            f.write('PLANE,' + fric_obj_name + '\n')

    def write_constraints_transform(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Transform Constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for trans_object in self.transform_objects:
            trans_obj = trans_object['Object']
            f.write('** ' + trans_obj.Label + '\n')
            if trans_obj.TransformType == "Rectangular":
                f.write('*TRANSFORM, NSET=Rect' + trans_obj.Name + ', TYPE=R\n')
                coords = FemMeshTools.get_rectangular_coords(trans_obj)
                f.write(coords + '\n')
            elif trans_obj.TransformType == "Cylindrical":
                f.write('*TRANSFORM, NSET=Cylin' + trans_obj.Name + ', TYPE=C\n')
                coords = FemMeshTools.get_cylindrical_coords(trans_obj)
                f.write(coords + '\n')

    def write_constraints_selfweight(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Self weight Constraint\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.selfweight_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            selwei_obj = femobj['Object']
            f.write('** ' + selwei_obj.Label + '\n')
            f.write('*DLOAD\n')
            f.write(self.ccx_eall + ',GRAV,9810,' + str(selwei_obj.Gravity_x) + ',' + str(selwei_obj.Gravity_y) + ',' + str(selwei_obj.Gravity_z) + '\n')
            f.write('\n')
        # grav (erdbeschleunigung) is equal for all elements
        # should be only one constraint
        # different element sets for different density are written in the material element sets already

    def write_constraints_force(self, f):
        # check shape type of reference shape and get node loads
        self.get_constraints_force_nodeloads()
        # write node loads to file
        f.write('\n***********************************************************\n')
        f.write('** Node loads Constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        f.write('*CLOAD\n')
        for femobj in self.force_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            f.write('** ' + femobj['Object'].Label + '\n')
            direction_vec = femobj['Object'].DirectionVector
            for ref_shape in femobj['NodeLoadTable']:
                f.write('** ' + ref_shape[0] + '\n')
                for n in sorted(ref_shape[1]):
                    node_load = ref_shape[1][n]
                    if (direction_vec.x != 0.0):
                        v1 = "{:.13E}".format(direction_vec.x * node_load)
                        f.write(str(n) + ',1,' + v1 + '\n')
                    if (direction_vec.y != 0.0):
                        v2 = "{:.13E}".format(direction_vec.y * node_load)
                        f.write(str(n) + ',2,' + v2 + '\n')
                    if (direction_vec.z != 0.0):
                        v3 = "{:.13E}".format(direction_vec.z * node_load)
                        f.write(str(n) + ',3,' + v3 + '\n')
                f.write('\n')
            f.write('\n')

    def write_constraints_pressure(self, f):
        # get the faces and face numbers
        self.get_constraints_pressure_faces()
        # write face loads to file
        f.write('\n***********************************************************\n')
        f.write('** Element + CalculiX face + load in [MPa]\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for femobj in self.pressure_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            prs_obj = femobj['Object']
            f.write('** ' + prs_obj.Label + '\n')
            rev = -1 if prs_obj.Reversed else 1
            f.write('*DLOAD\n')
            for ref_shape in femobj['PressureFaces']:
                f.write('** ' + ref_shape[0] + '\n')
                for face, fno in ref_shape[1]:
                    if fno > 0:  # solid mesh face
                        f.write("{},P{},{}\n".format(face, fno, rev * prs_obj.Pressure))
                    elif fno == 0:  # on shell mesh face: fno == 0 --> normal of element face == face normal
                        f.write("{},P,{}\n".format(face, rev * prs_obj.Pressure))
                    elif fno == -1:  # on shell mesh face: fno == -1 --> normal of element face opposite direction face normal
                        f.write("{},P,{}\n".format(face, -1 * rev * prs_obj.Pressure))

    def write_constraints_temperature(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Fixed temperature constraint applied\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for ftobj in self.temperature_objects:
            fixedtemp_obj = ftobj['Object']
            f.write('** ' + fixedtemp_obj.Label + '\n')
            NumberOfNodes = len(ftobj['Nodes'])
            if fixedtemp_obj.ConstraintType == "Temperature":
                f.write('*BOUNDARY\n')
                f.write('{},11,11,{}\n'.format(fixedtemp_obj.Name, fixedtemp_obj.Temperature))
                f.write('\n')
            elif fixedtemp_obj.ConstraintType == "CFlux":
                f.write('*CFLUX\n')
                f.write('{},11,{}\n'.format(fixedtemp_obj.Name, fixedtemp_obj.CFlux * 0.001 / NumberOfNodes))
                f.write('\n')

    def write_constraints_heatflux(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Heatflux constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        for hfobj in self.heatflux_objects:
            heatflux_obj = hfobj['Object']
            f.write('** ' + heatflux_obj.Label + '\n')
            if heatflux_obj.ConstraintType == "Convection":
                f.write('*FILM\n')
                for o, elem_tup in heatflux_obj.References:
                    for elem in elem_tup:
                        ho = o.Shape.getElement(elem)
                        if ho.ShapeType == 'Face':
                            v = self.mesh_object.FemMesh.getccxVolumesByFace(ho)
                            f.write("** Heat flux on face {}\n".format(elem))
                            for i in v:
                                f.write("{},F{},{},{}\n".format(i[0], i[1], heatflux_obj.AmbientTemp, heatflux_obj.FilmCoef * 0.001))  # SvdW add factor to force heatflux to units system of t/mm/s/K # OvG: Only write out the VolumeIDs linked to a particular face
            elif heatflux_obj.ConstraintType == "DFlux":
                f.write('*DFLUX\n')
                for o, elem_tup in heatflux_obj.References:
                    for elem in elem_tup:
                        ho = o.Shape.getElement(elem)
                        if ho.ShapeType == 'Face':
                            v = self.mesh_object.FemMesh.getccxVolumesByFace(ho)
                            f.write("** Heat flux on face {}\n".format(elem))
                            for i in v:
                                f.write("{},S{},{}\n".format(i[0], i[1], heatflux_obj.DFlux * 0.001))

    def write_constraints_fluidsection(self, f):
        f.write('\n***********************************************************\n')
        f.write('** FluidSection constraints\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        if os.path.exists(self.fluid_inout_nodes_file):
            inout_nodes_file = open(self.fluid_inout_nodes_file, "r")
            lines = inout_nodes_file.readlines()
            inout_nodes_file.close()
        else:
            FreeCAD.Console.PrintError("1DFlow inout nodes file not found: " + self.fluid_inout_nodes_file + '\n')
        # get nodes
        self.get_constraints_fluidsection_nodes()
        for femobj in self.fluidsection_objects:  # femobj --> dict, FreeCAD document object is femobj['Object']
            fluidsection_obj = femobj['Object']
            f.write('** ' + fluidsection_obj.Label + '\n')
            if fluidsection_obj.SectionType == 'Liquid':
                if fluidsection_obj.LiquidSectionType == 'PIPE INLET':
                    f.write('**Fluid Section Inlet \n')
                    if fluidsection_obj.InletPressureActive is True:
                        f.write('*BOUNDARY \n')
                        for n in femobj['Nodes']:
                            for line in lines:
                                b = line.split(',')
                                if int(b[0]) == n and b[3] == 'PIPE INLET\n':
                                    f.write(b[0] + ',2,2,' + str(fluidsection_obj.InletPressure) + '\n')  # degree of freedom 2 is for defining pressure
                    if fluidsection_obj.InletFlowRateActive is True:
                        f.write('*BOUNDARY,MASS FLOW \n')
                        for n in femobj['Nodes']:
                            for line in lines:
                                b = line.split(',')
                                if int(b[0]) == n and b[3] == 'PIPE INLET\n':
                                    f.write(b[1] + ',1,1,' + str(fluidsection_obj.InletFlowRate * 0.001) + '\n')  # degree of freedom 1 is for defining flow rate, factor applied to convert unit from kg/s to t/s
                elif fluidsection_obj.LiquidSectionType == 'PIPE OUTLET':
                    f.write('**Fluid Section Outlet \n')
                    if fluidsection_obj.OutletPressureActive is True:
                        f.write('*BOUNDARY \n')
                        for n in femobj['Nodes']:
                            for line in lines:
                                b = line.split(',')
                                if int(b[0]) == n and b[3] == 'PIPE OUTLET\n':
                                    f.write(b[0] + ',2,2,' + str(fluidsection_obj.OutletPressure) + '\n')  # degree of freedom 2 is for defining pressure
                    if fluidsection_obj.OutletFlowRateActive is True:
                        f.write('*BOUNDARY,MASS FLOW \n')
                        for n in femobj['Nodes']:
                            for line in lines:
                                b = line.split(',')
                                if int(b[0]) == n and b[3] == 'PIPE OUTLET\n':
                                    f.write(b[1] + ',1,1,' + str(fluidsection_obj.OutletFlowRate * 0.001) + '\n')  # degree of freedom 1 is for defining flow rate, factor applied to convert unit from kg/s to t/s

    def write_outputs_types(self, f):
        f.write('\n***********************************************************\n')
        f.write('** Outputs --> frd file\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        if self.beamsection_objects or self.shellthickness_objects or self.fluidsection_objects:
            if self.solver_obj.BeamShellResultOutput3D is False:
                f.write('*NODE FILE, OUTPUT=2d\n')
            else:
                f.write('*NODE FILE, OUTPUT=3d\n')
        else:
            f.write('*NODE FILE\n')
        if self.analysis_type == "thermomech":  # MPH write out nodal temperatures if thermomechanical
            if not self.fluidsection_objects:
                f.write('U, NT\n')
            else:
                f.write('MF, PS\n')
        else:
            f.write('U\n')
        if not self.fluidsection_objects:
            f.write('*EL FILE\n')
            if self.solver_obj.MaterialNonlinearity == 'nonlinear':
                f.write('S, E, PEEQ\n')
            else:
                f.write('S, E\n')
            # there is no need to write all integration point results as long as there is no reader for this
            # see https://forum.freecadweb.org/viewtopic.php?f=18&t=29060
            # f.write('** outputs --> dat file\n')
            # f.write('*NODE PRINT , NSET=' + self.ccx_nall + '\n')
            # f.write('U \n')
            # f.write('*EL PRINT , ELSET=' + self.ccx_eall + '\n')
            # f.write('S \n')

    def write_step_end(self, f):
        f.write('\n***********************************************************\n')
        f.write('** written by {} function\n'.format(sys._getframe().f_code.co_name))
        f.write('*END STEP \n')


    # self.ccx_elsets = [ {
    #                        'ccx_elset' : [e1, e2, e3, ... , en] or elements set name strings
    #                        'ccx_elset_name' : 'ccx_identifier_elset'
    #                        'mat_obj_name' : 'mat_obj.Name'
    #                        'ccx_mat_name' : 'mat_obj.Material['Name']'   !!! not unique !!!
    #                        'beamsection_obj' : 'beamsection_obj'         if exists
    #                        'fluidsection_obj' : 'fluidsection_obj'       if exists
    #                        'shellthickness_obj' : shellthickness_obj'    if exists
    #                        'beam_normal' : normal vector                 for beams only
    #                     },
    #                     {}, ... , {} ]

    # beam
    # TODO support multiple beamrotations
    # we do not need any more any data from the rotation document object, thus we do not need to save the rotation document object name in the else
    def get_ccx_elsets_single_mat_single_beam(self):
        mat_obj = self.material_objects[0]['Object']
        beamsec_obj = self.beamsection_objects[0]['Object']
        beamrot_data = self.beamrotation_objects[0]
        for i, beamdirection in enumerate(beamrot_data['FEMRotations1D']):
            elset_data = beamdirection['ids']  # ID's for this direction
            names = [{'short': 'M0'}, {'short': 'B0'}, {'short': beamrot_data['ShortName']}, {'short': 'D' + str(i)}]
            ccx_elset = {}
            ccx_elset['ccx_elset'] = elset_data
            ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
            ccx_elset['mat_obj_name'] = mat_obj.Name
            ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
            ccx_elset['beamsection_obj'] = beamsec_obj
            ccx_elset['beam_normal'] = beamdirection['normal']  # normal for this direction
            self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_single_mat_multiple_beam(self):
        mat_obj = self.material_objects[0]['Object']
        beamrot_data = self.beamrotation_objects[0]
        for beamsec_data in self.beamsection_objects:
            beamsec_obj = beamsec_data['Object']
            beamsec_ids = set(beamsec_data['FEMElements'])
            for i, beamdirection in enumerate(beamrot_data['FEMRotations1D']):
                beamdir_ids = set(beamdirection['ids'])
                elset_data = list(sorted(beamsec_ids.intersection(beamdir_ids)))  # empty intersection sets possible
                if elset_data:
                    names = [{'short': 'M0'}, {'short': beamsec_data['ShortName']}, {'short': beamrot_data['ShortName']}, {'short': 'D' + str(i)}]
                    ccx_elset = {}
                    ccx_elset['ccx_elset'] = elset_data
                    ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
                    ccx_elset['mat_obj_name'] = mat_obj.Name
                    ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
                    ccx_elset['beamsection_obj'] = beamsec_obj
                    ccx_elset['beam_normal'] = beamdirection['normal']  # normal for this direction
                    self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_multiple_mat_single_beam(self):
        beamsec_obj = self.beamsection_objects[0]['Object']
        beamrot_data = self.beamrotation_objects[0]
        for mat_data in self.material_objects:
            mat_obj = mat_data['Object']
            mat_ids = set(mat_data['FEMElements'])
            for i, beamdirection in enumerate(beamrot_data['FEMRotations1D']):
                beamdir_ids = set(beamdirection['ids'])
                elset_data = list(sorted(mat_ids.intersection(beamdir_ids)))
                if elset_data:
                    names = [{'short': mat_data['ShortName']}, {'short': 'B0'}, {'short': beamrot_data['ShortName']}, {'short': 'D' + str(i)}]
                    ccx_elset = {}
                    ccx_elset['ccx_elset'] = elset_data
                    ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
                    ccx_elset['mat_obj_name'] = mat_obj.Name
                    ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
                    ccx_elset['beamsection_obj'] = beamsec_obj
                    ccx_elset['beam_normal'] = beamdirection['normal']  # normal for this direction
                    self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_multiple_mat_multiple_beam(self):
        beamrot_data = self.beamrotation_objects[0]
        for beamsec_data in self.beamsection_objects:
            beamsec_obj = beamsec_data['Object']
            beamsec_ids = set(beamsec_data['FEMElements'])
            for mat_data in self.material_objects:
                mat_obj = mat_data['Object']
                mat_ids = set(mat_data['FEMElements'])
                for i, beamdirection in enumerate(beamrot_data['FEMRotations1D']):
                    beamdir_ids = set(beamdirection['ids'])
                    elset_data = list(sorted(beamsec_ids.intersection(mat_ids).intersection(beamdir_ids)))  # empty intersection sets possible
                    if elset_data:
                        names = [{'short': mat_data['ShortName']}, {'short': beamsec_data['ShortName']}, {'short': beamrot_data['ShortName']}, {'short': 'D' + str(i)}]
                        ccx_elset = {}
                        ccx_elset['ccx_elset'] = elset_data
                        ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
                        ccx_elset['mat_obj_name'] = mat_obj.Name
                        ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
                        ccx_elset['beamsection_obj'] = beamsec_obj
                        ccx_elset['beam_normal'] = beamdirection['normal']  # normal for this direction
                        self.ccx_elsets.append(ccx_elset)

    # fluid
    def get_ccx_elsets_single_mat_single_fluid(self):
        mat_obj = self.material_objects[0]['Object']
        fluidsec_obj = self.fluidsection_objects[0]['Object']
        elset_data = self.ccx_eedges
        names = [{'short': 'M0'}, {'short': 'F0'}]
        ccx_elset = {}
        ccx_elset['ccx_elset'] = elset_data
        ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
        ccx_elset['mat_obj_name'] = mat_obj.Name
        ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
        ccx_elset['fluidsection_obj'] = fluidsec_obj
        self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_single_mat_multiple_fluid(self):
        mat_obj = self.material_objects[0]['Object']
        for fluidsec_data in self.fluidsection_objects:
            fluidsec_obj = fluidsec_data['Object']
            elset_data = fluidsec_data['FEMElements']
            names = [{'short': 'M0'}, {'short': fluidsec_data['ShortName']}]
            ccx_elset = {}
            ccx_elset['ccx_elset'] = elset_data
            ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
            ccx_elset['mat_obj_name'] = mat_obj.Name
            ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
            ccx_elset['fluidsection_obj'] = fluidsec_obj
            self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_multiple_mat_single_fluid(self):
        fluidsec_obj = self.fluidsection_objects[0]['Object']
        for mat_data in self.material_objects:
            mat_obj = mat_data['Object']
            elset_data = mat_data['FEMElements']
            names = [{'short': mat_data['ShortName']}, {'short': 'F0'}]
            ccx_elset = {}
            ccx_elset['ccx_elset'] = elset_data
            ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
            ccx_elset['mat_obj_name'] = mat_obj.Name
            ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
            ccx_elset['fluidsection_obj'] = fluidsec_obj
            self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_multiple_mat_multiple_fluid(self):
        for fluidsec_data in self.fluidsection_objects:
            fluidsec_obj = fluidsec_data['Object']
            for mat_data in self.material_objects:
                mat_obj = mat_data['Object']
                fluidsec_ids = set(fluidsec_data['FEMElements'])
                mat_ids = set(mat_data['FEMElements'])
                elset_data = list(sorted(fluidsec_ids.intersection(mat_ids)))  # empty intersection sets possible
                if elset_data:
                    names = [{'short': mat_data['ShortName']}, {'short': fluidsec_data['ShortName']}]
                    ccx_elset = {}
                    ccx_elset['ccx_elset'] = elset_data
                    ccx_elset['ccx_elset_name'] = get_ccx_elset_name_short(names)
                    ccx_elset['mat_obj_name'] = mat_obj.Name
                    ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
                    ccx_elset['fluidsection_obj'] = fluidsec_obj
                    self.ccx_elsets.append(ccx_elset)

    # shell
    def get_ccx_elsets_single_mat_single_shell(self):
        mat_obj = self.material_objects[0]['Object']
        shellth_obj = self.shellthickness_objects[0]['Object']
        elset_data = self.ccx_efaces
        names = [{'long': mat_obj.Name, 'short': 'M0'}, {'long': shellth_obj.Name, 'short': 'S0'}]
        ccx_elset = {}
        ccx_elset['ccx_elset'] = elset_data
        ccx_elset['ccx_elset_name'] = get_ccx_elset_name_standard(names)
        ccx_elset['mat_obj_name'] = mat_obj.Name
        ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
        ccx_elset['shellthickness_obj'] = shellth_obj
        self.ccx_elsets.append(ccx_elset)

 
    # solid
    def get_ccx_elsets_single_mat_solid(self):
        mat_obj = self.material_objects[0]['Object']
        elset_data = self.ccx_evolumes
        names = [{'long': mat_obj.Name, 'short': 'M0'}, {'long': 'Solid', 'short': 'Solid'}]
        ccx_elset = {}
        ccx_elset['ccx_elset'] = elset_data
        ccx_elset['ccx_elset_name'] = get_ccx_elset_name_standard(names)
        ccx_elset['mat_obj_name'] = mat_obj.Name
        ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
        self.ccx_elsets.append(ccx_elset)

    def get_ccx_elsets_multiple_mat_solid(self):
        for mat_data in self.material_objects:
            mat_obj = mat_data['Object']
            elset_data = mat_data['FEMElements']
            names = [{'long': mat_obj.Name, 'short': mat_data['ShortName']}, {'long': 'Solid', 'short': 'Solid'}]
            ccx_elset = {}
            ccx_elset['ccx_elset'] = elset_data
            ccx_elset['ccx_elset_name'] = get_ccx_elset_name_standard(names)
            ccx_elset['mat_obj_name'] = mat_obj.Name
            ccx_elset['ccx_mat_name'] = mat_obj.Material['Name']
            self.ccx_elsets.append(ccx_elset)



##  @}
