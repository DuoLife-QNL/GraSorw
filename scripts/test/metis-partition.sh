program_metis_parition=/home/hongzheng/Codes/CLionProjects/IOE-SORW/metis-5.1.0/build/Linux-x86_64/programs/gpmetis

program_tometis=/home/hongzheng/Codes/CLionProjects/IOE-SORW/cmake-build-idmg-release/MetisPreProcess

LJ=/home/hongzheng/data/soc-LiveJournal1.txt.unDir
TW=/home/hongzheng/data/twitter_rv.net.unDir
UK=/home/hongzheng/data/uk_edge.txt.unDir
FR=/home/hongzheng/data/com-friendster.txt.unDir

#${program_tometis} -f ${LJ}
#${program_tometis} -f ${UK}
#${program_tometis} -f ${FR}
#
#${program_metis_parition} ${LJ}.toMetis 17
#${program_metis_parition} ${UK}.toMetis 25
${program_metis_parition} ${FR}.toMetis 27

