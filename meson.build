project(
    'class-group-tabulation',
    ['c', 'cpp'],
    default_options: ['c_std=gnu99', 'warning_level=2', 'optimization=3'],
    meson_version: '>=1.1.0',
)

cc = meson.get_compiler('c')
gmp_dep = dependency('gmp', required: true)
mpi_dep = dependency('mpi', language: 'c', required: true)
m_dep = cc.find_library('m', required: true)

with_pari = get_option('with_pari')
pari_deps = []
if with_pari
    pari_deps += [cc.find_library('pari', required: true)]
endif

vendored_target_options = ['warning_level=0']

build_root = meson.current_build_dir()
source_root = meson.current_source_dir()
staged_include_root = join_paths(build_root, 'staged-include')
generated_include_root = join_paths(build_root, 'generated-include')

foreach dir_path : [
    staged_include_root,
    generated_include_root,
    join_paths(generated_include_root, 'liboptarith'),
]
    run_command('mkdir', '-p', dir_path, check: true)
endforeach

foreach staged_link : [
    [
        join_paths(source_root, 'libs/optarith'),
        join_paths(staged_include_root, 'liboptarith'),
    ],
    [
        join_paths(source_root, 'libs/qform'),
        join_paths(staged_include_root, 'libqform'),
    ],
    [
        join_paths(source_root, 'libs/clgrp'),
        join_paths(staged_include_root, 'libclgrp'),
    ],
    [
        join_paths(build_root, 'sqrtmodp_list.h'),
        join_paths(generated_include_root, 'liboptarith', 'sqrtmodp_list.h'),
    ],
]
    run_command('ln', '-sfn', staged_link[0], staged_link[1], check: true)
endforeach

gen_sqrtmodp = executable(
    'gen_sqrtmodp',
    ['libs/optarith/code_gen/gen_sqrtmodp.cc', 'libs/optarith/primes.c'],
    include_directories: include_directories('libs/optarith'),
    dependencies: [gmp_dep],
    c_args: ['-w'],
    cpp_args: ['-w'],
    override_options: vendored_target_options,
    native: true,
    build_by_default: false,
)

sqrtmodp_generated = custom_target(
    'sqrtmodp_generated',
    output: ['sqrtmodp_list.c', 'sqrtmodp_list.h'],
    command: [
        gen_sqrtmodp,
        get_option('sqrtmodp_maxp').to_string(),
        '@OUTPUT0@',
        '@OUTPUT1@',
    ],
)

sqrtmodp_c = sqrtmodp_generated[0]
sqrtmodp_h = sqrtmodp_generated[1]

vendored_include_args = [
    '-isystem',
    generated_include_root,
    '-isystem',
    staged_include_root,
]
vendored_compile_args = vendored_include_args + ['-w']
app_compile_args = vendored_include_args + [
    '-I' + join_paths(source_root, 'src'),
]

optarith_sources = files(
    'libs/optarith/closest_23.c',
    'libs/optarith/gcd/gcd_binary_l2r.c',
    'libs/optarith/gcd/gcd_brent.c',
    'libs/optarith/gcd/gcd_common.c',
    'libs/optarith/gcd/gcd_divrem.c',
    'libs/optarith/gcd/gcd_flint.c',
    'libs/optarith/gcd/gcd_lehmer.c',
    'libs/optarith/gcd/gcd_mpz128.c',
    'libs/optarith/gcd/gcd_shallit.c',
    'libs/optarith/gcd/gcd_smallq_case.c',
    'libs/optarith/gcd/gcd_smallq_loop.c',
    'libs/optarith/gcd/gcd_stein.c',
    'libs/optarith/group.c',
    'libs/optarith/group_pow.c',
    'libs/optarith/lehmer_table.c',
    'libs/optarith/math32.c',
    'libs/optarith/math64.c',
    'libs/optarith/math_mpz.c',
    'libs/optarith/mpz_xgcd.c',
    'libs/optarith/primes.c',
    'libs/optarith/primorial.c',
    'libs/optarith/s128_t.c',
    'libs/optarith/square_free.c',
    'libs/optarith/u128_t.c',
) + [sqrtmodp_c]

optarith = static_library(
    'optarith',
    optarith_sources + [sqrtmodp_h],
    c_args: vendored_compile_args,
    dependencies: [gmp_dep, m_dep],
    override_options: vendored_target_options,
)

optarith_dep = declare_dependency(
    link_with: optarith,
    dependencies: [gmp_dep, m_dep],
)

qform = static_library(
    'qform',
    files(
        'libs/qform/costs_128.c',
        'libs/qform/costs_64.c',
        'libs/qform/dbreps/mpz_pow_reps.c',
        'libs/qform/dbreps/s128_pow_reps.c',
        'libs/qform/dbreps/s64_pow_reps.c',
        'libs/qform/gen_qform.c',
        'libs/qform/mpz_qform.c',
        'libs/qform/qform_group.c',
        'libs/qform/s128_qform.c',
        'libs/qform/s64_qform.c',
    ) + [sqrtmodp_h],
    c_args: vendored_compile_args,
    dependencies: [gmp_dep],
    override_options: vendored_target_options,
)

qform_dep = declare_dependency(
    link_with: qform,
    dependencies: [optarith_dep, gmp_dep, m_dep],
)

clgrp = static_library(
    'clgrp',
    files(
        'libs/clgrp/clgrp.c',
        'libs/clgrp/functions.c',
        'libs/clgrp/sieve.c',
        'libs/clgrp/verify.c',
    ) + [sqrtmodp_h],
    c_args: vendored_compile_args + ['-DKEEP_FILES'],
    dependencies: [gmp_dep],
    override_options: vendored_target_options,
)

clgrp_dep = declare_dependency(
    link_with: clgrp,
    dependencies: [qform_dep, optarith_dep, gmp_dep, m_dep],
)

clgrp_ell = executable(
    'clgrp_ell',
    files('src/clgrp_ell.c', 'src/clgrp_ell_main.c'),
    c_args: app_compile_args + (with_pari ? ['-DWITH_PARI'] : []),
    dependencies: [clgrp_dep, mpi_dep] + pari_deps,
)

ell_count = executable(
    'ell_count',
    'src/ell_count_main.c',
    c_args: app_compile_args,
    dependencies: [mpi_dep],
)

alias_target('apps', [clgrp_ell, ell_count])

summary(
    {
        'with_pari': with_pari,
        'sqrtmodp_maxp': get_option('sqrtmodp_maxp'),
    },
)
