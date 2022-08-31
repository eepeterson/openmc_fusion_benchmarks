import openmc.deplete
import openmc

model = openmc.Model.from_xml()
materials = model.geometry.get_all_materials()
nuclides = set()
for material in materials.values():
    nuclides |= set(material.get_nuclides())
nuclides.remove('C0')
nuclides |= {'C12', 'C13'}

chain = openmc.deplete.Chain.from_xml('/home/romano/openmc-data/depletion/chain_endfb71_pwr.xml')
reduced_chain = chain.reduce(nuclides, 4)
print(len(reduced_chain))
reduced_chain.export_to_xml('chain_reduced.xml')
