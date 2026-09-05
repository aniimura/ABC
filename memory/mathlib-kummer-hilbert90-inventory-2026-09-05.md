---
name: mathlib-kummer-hilbert90-inventory-2026-09-05
description: Kummer は 1 本の巡回拡大のみ。Hilbert 90 は有限次必須。X^n-C a の既約判定は奇数 n 限定
metadata:
  type: reference
---

pin mathlib `db127794` + ABC3 の実測(2026-09-05)。pGC 経路 C(`ResearchPaper/pgc-goal.md`)向け。

## ★落とし穴 2 つ

1. **Hilbert 90 は `[FiniteDimensional K L]` 必須**。
   `groupCohomology.H1ofAutOnUnitsUnique` /
   `isMulCoboundary₁_of_isMulCocycle₁_of_aut_to_units` /
   `exists_div_of_norm_eq_one`(`RepresentationTheory/Homological/GroupCohomology/Hilbert90.lean`)
   ——4 つの variable ブロック全部に付いており、証明の核が `Gal(L/K)` 上の**有限和**。
   ファイルに「Develop Galois cohomology to extend Noether's result to infinite Galois
   extensions」と TODO 明記。**`Γ_K` や `K^ur/K` には直接使えない。**
2. **`X^n − C a` の既約判定は奇数 `n` 限定**。
   `X_pow_sub_C_irreducible_iff_of_odd` / `..._of_prime_pow (hp' : p ≠ 2)` のみで、
   ファイルに `TODO: criteria for even n`。**`p = 2` が素通しできない。**

## Kummer 理論(`Mathlib.FieldTheory.KummerExtension`)

**ある**: `autEquivZmod (H) (L) (hζ) : Gal(L/K) ≃* Multiplicative (ZMod n)`、
`autEquivRootsOfUnity`、`isCyclic_tfae`、`exists_root_adjoin_eq_top_of_isCyclic`。
★`μ_n ⊆ K` は型クラスでなく明示引数 `hζ : (primitiveRoots n K).Nonempty`。
`n` は拡大次数に固定。引数順が不規則(`H` が先、`hζ` が後)。

**無い**: **`K^×/(K^×)^n ↔ 指数 n アーベル拡大`の大域双対**。
`KummerExtension.lean` を全文直読して 0 件。`Kˣ ⧸ (powMonoidHom n).range` を
Galois 群に結ぶ宣言は存在しない。`KummerDedekind` は素イデアル分解で無関係。

## 数え上げの部品

- ★**`IsCyclic.index_powMonoidHom_range (d) : (powMonoidHom d).range.index = (Nat.card G).gcd d`**
  (`GroupTheory/SpecificGroups/Cyclic.lean:635`)——**pGC 経路 C-q の `gcd(n, q−1)` そのもの**。
  ABC3 の `isCyclic_primeToPTorsion` + `card_primeToPTorsion` と合わせれば即出る。
- `Subgroup.quotient_finite_of_isOpen` / `index_map_equiv` / `index_ker` / `index_prod` / `index_pi`
- `QuotientGroup.mulEquivPiModRangePowMonoidHom`(積の `A/A^n`)
- `CommGroup.card_monoidHom_of_hasEnoughRootsOfUnity`
  ——★行き先が `Mˣ` で `[HasEnoughRootsOfUnity M (exponent G)]` 必須。`M = ZMod n` には当たらない。
- **`ContinuousMonoidHom` に濃度・`ker`・`range` の宣言は 1 つも無い**(98 宣言を全列挙)。
  `Hom_cont(Γ, ℤ/n)` の数え上げは自作(3 行:`ker` が開 ⇒ 商が有限 ⇒ `Abelianization.lift`)。

## Hensel(ABC3 の自前が主力)

- `ABC3.Found.exists_root_of_residue_root`(`Found/HenselianSplits.lean:46`、単根の持ち上げ)
- ★**`ABC3.Found.teichmuller : (IsLocalRing.ResidueField A)ˣ →* Aˣ`**(`Found/Teichmuller.lean`)
  ——「`q−1` 乗根がすべて持ち上がる」の実体。`n ∣ q−1` の `n` 乗根はここから直接出る。
- `henselianLocalRing_carrierIntegers` / `henselianLocalRing_adjoinIntegers`(instance 済み)
- mathlib 側は `HenselianLocalRing`(★`IsHenselianLocalRing` ではない)のみ。
  **「`p ∤ n` なら局所環の単数は `n` 乗数」は mathlib にも ABC3 にも無い。**
  `RootableBy` はあるが Henselian 局所環の単数群への instance は無い。

## ★`𝒪_{K^ur}` は ABC3 に無い

`unramifiedClosure K : IntermediateField` はあるが、整数環は**有限部分拡大ごとの
`adjoinIntegers K x` しか無い**。`HenselianLocalRing` もそちらにしか付いていない。
→ 「`𝒪_{K^ur}^×` は n 可除」は**有限段で述べて合併を取る**のが安い。

## ★重複の発見

`ABC3.Found.PGC.map_powRange_of_mulEquiv`(`Found/PGC/UnitsGroupInvariants.lean:78`)は
mathlib の `MulEquiv.map_range_powMonoidHom`(`Algebra/Group/Subgroup/Ker.lean:602`)と
同じ主張。`index_powRange_of_mulEquiv` も `Subgroup.index_map_equiv` から 1 行で出る。

## Herbrand 商

mathlib・ABC3 とも **0 件**。ただし経路 C は Herbrand 商を要求しない
(必要なのは `|A/A^n|` の同型不変性だけで、上の移送補題で足りる)。
