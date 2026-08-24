/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import ABC3.Found.Divisor.HeightOneDVR

/-!
# 代数的 Hartogs —— 正規 Noether 整域は高さ 1 の局所化の交わり(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV,

## ★★なぜ要るか —— `Example 6.1` に残る**唯一の壁**の片割れ

`Example 6.1` は model Frobenioid `C_{V,K̄,D_K}` を作るところまで済んでいる
(`ex61-model`)。残るのは原文が "(we observe that)" で畳んだ 1 行

```
𝒪^×(A) = 𝒪^▷(A) = k_L^×
```

だけである(節点 `ex61-units` / `global-units`)。★模型 Frobenioid では
`𝒪^×(A) = {u ∈ B(L) | div_{D_L}(u) = 0}` であり、`B(L)` の定義から
`D_L` の外では `ord = 0` なので、これは
**「すべての余次元 1 の点で `ord_v(u) = 0`」**と同値である。

★そこから `u ∈ k_L^×` を出すのに 2 本要る:

| 要るもの | 状態 |
|---|---|
| **代数的 Hartogs**(本ファイル) | ★mathlib に無い(`IsKrullDomain` 不在) |
| `Γ(X, 𝒪_X)` が `k` 上有限な体(proper) | ★★**mathlib にあった**(下記) |

## ★★★★★在庫の訂正(2026-08-25、実測)

依存グラフは `proper-global-sections` を `absent` と記録していたが、
★★**これは誤りである**。`Mathlib/AlgebraicGeometry/Morphisms/Proper.lean` に

* `AlgebraicGeometry.isField_of_universallyClosed`
  —— `X` が整で `Spec K` 上 universally closed なら `Γ(X, ⊤)` は**体**
* `AlgebraicGeometry.finite_appTop_of_universallyClosed`
  —— さらに locally of finite type なら `Γ(X, ⊤)` は `K` 上**有限**

があり、`IsProper = IsSeparated ∧ UniversallyClosed ∧ LocallyOfFiniteType` なので
**proper ならそのまま当たる**。★したがって `Example 6.1` に残る在庫不足は
**本ファイルの代数的 Hartogs 1 本だけ**である。

## ★中身(紙の上)

`R` が Noether 整閉整域、`K = Frac R` とする。`x ∈ K` について

```
x ∈ R  ⟺  ∀ v : 高さ 1 の素, x ∈ R_v
```

★「⟸」が本体。分母イデアル `I = {r ∈ R | r·x ∈ R}` を取ると、
`I ≠ R` なら `I` の随伴素 `p` が取れる。`R` が正規なら
**`p` の高さは 1**(Serre の判定 `R1 + S2` の `S2` の側)なので、
仮定から `x ∈ R_p` となり `I ⊄ p` に矛盾する。

★★mathlib の在庫で近いのは `MaximalSpectrum.iInf_localization_eq_bot`
(**極大**イデアルにわたる交わり)だけで、**高さ 1 版ではない**。
-/

namespace ABC3.Skeleton.Hartogs

open ABC3.Meta ABC3.Found.Divisor

universe u

variable {R : Type u} [CommRing R] [IsDomain R] [IsNoetherianRing R] [IsIntegrallyClosed R]

/-- ★★★★★**代数的 Hartogs** ——
正規 Noether 整域の元は、**すべての高さ 1 の素での局所化に入ること**で特徴づけられる。

★`x ∈ R_v` を「`v` での位数が非負」ではなく
**`∃ r s, x = r/s ∧ s ∉ v`** の形で書く(局所化の API 摩擦を避ける)。 -/
theorem mem_of_forall_heightOne (x : FractionRing R)
    (_h : ∀ v : HeightOnePrime R, ∃ r s : R, s ∉ v.asIdeal ∧
      x * algebraMap R (FractionRing R) s = algebraMap R (FractionRing R) r) :
    ∃ a : R, x = algebraMap R (FractionRing R) a := by
  sorry

/-- ★★**逆向き**(自明) —— `R` の元はどの局所化にも入る。 -/
theorem forall_heightOne_of_mem (a : R) (v : HeightOnePrime R) :
    ∃ r s : R, s ∉ v.asIdeal ∧
      algebraMap R (FractionRing R) a * algebraMap R (FractionRing R) s
        = algebraMap R (FractionRing R) r :=
  ⟨a, 1, fun h => v.isPrime.ne_top (Ideal.eq_top_of_isUnit_mem _ h isUnit_one), by simp⟩

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def mem_of_forall_heightOne.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 代数的 Hartogs(正規 Noether 整域は高さ 1 の局所化の交わり)",
    sectionId := "frdi-example-6-1" }

/-- ★★★**この節点が §6 の律速の 1 本である**。

★訂正(2026-08-25): 相方の `proper-global-sections` は
**mathlib に在庫があった**(`isField_of_universallyClosed` /
`finite_appTop_of_universallyClosed`)。したがって
`Example 6.1` の条なし `.src` を止めているのは**本節点だけ**である。 -/
def mem_of_forall_heightOne.needs : List ProofObligation :=
  [ .citation "[mathlib]" "MaximalSpectrum.iInf_localization_eq_bot(極大イデアルにわたる交わり)"
      (.inMathlib "MaximalSpectrum.iInf_localization_eq_bot") 110,
    .citation "[mathlib]" "IsKrullDomain / 高さ 1 の素にわたる交わり"
      (.absent "Mathlib/RingTheory/ に KrullDomain は無い(KrullDimension のみ)。LocalProperties/IntegrallyClosed.lean にも高さ 1 版は無い(2026-08-25 に ls で実測)") 110,
    .citation "[mathlib]" "AlgebraicGeometry.isField_of_universallyClosed(proper なら Γ(X,⊤) は体)"
      (.inMathlib "AlgebraicGeometry.isField_of_universallyClosed") 110,
    .citation "[mathlib]" "AlgebraicGeometry.finite_appTop_of_universallyClosed(Γ(X,⊤) は k 上有限)"
      (.inMathlib "AlgebraicGeometry.finite_appTop_of_universallyClosed") 110,
    .citation "[ABC3]" "isDiscreteValuationRing_localization_atPrime_of_minimal(高さ 1 で DVR)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isDiscreteValuationRing_localization_atPrime_of_minimal") 110,
    .derivation
      "分母イデアル I = {r | r·x ∈ R} を取り、I ≠ R なら随伴素 p が取れる。正規性(S2)から ht p = 1、仮定より x ∈ R_p、よって I ⊄ p に矛盾" 110,
    .implicitStep
      "★原文は「(we observe that)」の 1 語で畳んでいる。畳まれているのは代数的 Hartogs と proper 大域切断の 2 本" 110 ]

def forall_heightOne_of_mem.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 代数的 Hartogs の自明な向き",
    sectionId := "frdi-example-6-1" }

/-- ★自明な向き —— 依存は無い。 -/
def forall_heightOne_of_mem.needs : List ProofObligation := []

end ABC3.Skeleton.Hartogs
