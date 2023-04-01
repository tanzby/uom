#!/bin/bash

# Fetch all third-party libs.
bazel fetch "//uom/..."

# Append bazel_compdb plugin config to The WORKSPACE file.
readonly WORKSPACE="$(bazel info workspace)"
readonly WORKSPACE_FILE="${WORKSPACE}/WORKSPACE"
readonly BAZEL_COMPDB='
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
git_repository(
    name = "com_grail_bazel_compdb",
    remote = "https://github.com/grailbio/bazel-compilation-database.git",
    tag = "0.5.0",
)
'
WORKSPACE_FILE_CONTENT=`cat ${WORKSPACE_FILE}`
echo "${BAZEL_COMPDB}" >> "${WORKSPACE_FILE}"
 
# Create the BUILD file to generate compilation_database.json.
readonly BUILD_FILE_CONTENT='
load("@com_grail_bazel_compdb//:aspects.bzl", "compilation_database")
compilation_database(
    name = "gen",
    targets = [
        "//uom/env_examples:sophus",
        "//uom/env_examples:publish_rviz_marker",
        "//uom/env_examples:rosbag_reader",
    ],
)
'
mkdir "${WORKSPACE}/gen_compile_commands"
echo "${BUILD_FILE_CONTENT}" > "${WORKSPACE}/gen_compile_commands/BUILD"

# Command to generate the compilation database file.
bazel build "//gen_compile_commands:gen"

# Remove the BUILD file.
rm -rf "${WORKSPACE}/gen_compile_commands"

# Restore the WORKSPACE file.
echo "$WORKSPACE_FILE_CONTENT" > "${WORKSPACE_FILE}"

# Location of the compilation database file.
readonly OUTFILE="$(bazel info bazel-bin)/gen_compile_commands/compile_commands.json"

# Command to replace the marker for exec_root in the file.
sed -i "s@__EXEC_ROOT__@${WORKSPACE}@g" "${OUTFILE}"
# Command to replace -isystem to -I.
sed -i "s@-isystem external@-I external@g" "${OUTFILE}"

# Link the external fold to workspace.
readonly OUTPUT_BASE="$(bazel info output_base)"
ln -f -s "${OUTPUT_BASE}/external" "${WORKSPACE}"

# Copy to workspace.
rm -f "${WORKSPACE}/compile_commands.json"
cp ${OUTFILE} ${WORKSPACE}

# The compilation database is now ready to use at this location.
echo "Compilation Database: ${OUTFILE}"
