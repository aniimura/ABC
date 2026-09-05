---
name: mathlib-cohomology-inventory-2026-09-05
description: pin mathlib(db127794)の Galois コホモロジー在庫。連続版は入ったがカップ積も Brauer 群の乗算も無い
metadata:
  type: reference
---

pin されている mathlib rev `db127794` の実測(2026-09-05、`.cache/mathlib-index.txt` 232,023 宣言)。

**ある**

- `continuousCohomology`(`Mathlib.Algebra.Category.ContinuousCohomology.Basic`、15 宣言)
  ——`[IsTopologicalGroup G]` のみ要求。副有限性は不要なので Γ_K(krullTopology)に
  型としては当てはまる。★**ただし `H⁰ ≅ 不変元` 以外に定理が無い**。離散版との一致も
  長完全列も TODO。係数は `Rep k G` でなく `Action (TopModuleCat R) G`。
- `cyclotomicCharacter`(`Mathlib.NumberTheory.Cyclotomic.CyclotomicCharacter`)
  ——`spec`・`modularCyclotomicCharacter`・★**`cyclotomicCharacter.continuous`**(Krull 位相に対する連続性)まで揃う。
- `MonoidHom.transfer`(`Mathlib.GroupTheory.Transfer`)——`[H.FiniteIndex]` 必須、
  値域は `CommGroup`。`Verlagerung` では引けない。
- `ProfiniteGrp` / `InfiniteGalois.profiniteGalGrp`(副有限群としての Γ_K の土台)

**無い**

- **局所 Tate 双対性**。`tateDuality`/`poitou` は 0 件。★**カップ積が 0 件**なので、
  対 `Hⁱ × H²⁻ⁱ → H²` そのものが未定義——「述べることすらできない」段階。
- **`Br(K) ≅ ℚ/ℤ` の不変写像**。`Algebra/BrauerGroup/` は `Defs.lean` の 7 宣言のみで、
  ★`BrauerGroup K` には**群構造すら入っていない**(Setoid 商の型)。crossed product も 0 件。
- **`ℤ_p(1)` を Galois 加群として表す語彙**。`tateTwist`/`tateModule` は mathlib に 0 件。
  `cyclotomicCharacter` は自ファイル以外のどの statement にも現れない。
- `groupCohomology` は今も**離散**(`Basic.lean` の TODO に `Profinite cohomology.` と明記)。
  `tateCohomology` は `[Fintype G]` 必須。

**How to apply:** pGC Proposition 1.1 の「② 局所 Tate 双対性が未構築」経路は、
カップ積が無い時点で当分閉じない。距離が縮まったのは連続コホモロジーの**定義**だけで、
局所体特有の算術内容(不変写像)はまったく縮まっていない。
[[pgc-prop12-reciprocity-gap]] と同じ結論。次に測るなら
「`BrauerGroup` に乗算が入ったか」が最短の指標。
