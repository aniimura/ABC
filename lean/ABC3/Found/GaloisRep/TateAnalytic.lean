import ABC3.Found.GaloisRep.BasisFree
import Mathlib.Topology.Algebra.InfiniteSum.Nonarchimedean
import Mathlib.Topology.Algebra.Valued.ValuationTopology

/-!
# Galois (G6) 第 211 ブロック —— **★★★★★Tate 級数の解析の背骨は mathlib に在った**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★(G6) の在庫を測った(2026-08-20)

(G5)(G7)(G8) はいずれも `TateCurveData` を台に持つ——★**(G6) が Galois 側の残り 4 欄すべての
律速**である。そこで Tate 一意化 `E(K) ≅ Kˣ/q^ℤ` の在庫を測った。

| 段 | 内容 | mathlib 在庫 |
|---|---|---|
| 収束の判定 | 非アルキメデス完備群で `Summable f ↔ f → 0`(cofinite) | ★★`NonarchimedeanAddGroup.summable_iff_tendsto_cofinite_zero` ✓ |
| 級数の積 | `(Σf)(Σg) = Σ f_i g_j` | ★★`HasSum.mul_of_nonarchimedean` ✓ |
| 進位相 | 完備 DVR の `𝔪` 進位相は非アルキメデス | ★`AdicTopology` に `instance : NonarchimedeanRing R` ✓ |
| **付値位相** | `Valued K Γ₀` から非アルキメデス性 | ★**無い**(2026-08-20 実測、`infer_instance` が失敗) |
| Tate 級数そのもの | `a₄(q), a₆(q)`、`X(u,q), Y(u,q)` | ★**0 件**(`Tate curve` で全文検索して 0 件) |
| 群同型の検証 | `Kˣ/q^ℤ ≅ E_q(K)` | ★**0 件** |

★★★**当初「Tate 曲線は無い」とだけ測っていたが、解析の背骨——非アルキメデス総和法——は
在った**。★★欠けているのは Tate 級数そのものと群同型の検証であり、
「解析の基礎から積む」必要は無い。

## ★★本ブロックで埋める穴

`Valued K Γ₀` から `NonarchimedeanAddGroup K` が**出ない**(instance が無い)。
★付値の球 `{x | v x < γ}` は開加法部分群なので、これは 10 行で埋まる。
★★埋めれば `summable_iff_tendsto_cofinite_zero` が付値体でそのまま使える。

★★★グローバル instance にはしない——`Valued` を持つすべての環に当たると
探索が重くなるため、`haveI` で局所的に入れる形にした。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `nonarchimedeanAddGroup_of_valued` | ★★★**付値位相は非アルキメデス** |
| `summable_of_valued_tendsto_zero` | ★★★★**完備な付値体では項が 0 に行けば総和可能** |
-/

namespace ABC3.Found.GaloisRep

open Filter Topology

/-- ★★★**付値位相をもつ環は非アルキメデスである**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★付値の球 `{x | v x < γ}` が開加法部分群であることに尽きる。
★★グローバル instance にはしない(探索が重くなる)——`haveI` で使う。 -/
theorem nonarchimedeanAddGroup_of_valued {K Γ₀ : Type}
    [Ring K] [LinearOrderedCommGroupWithZero Γ₀] [Valued K Γ₀] :
    NonarchimedeanAddGroup K where
  is_nonarchimedean U hU := by
    obtain ⟨γ, hγ⟩ := Valued.mem_nhds_zero.1 hU
    let H : AddSubgroup K :=
      { carrier := {x : K | Valued.v.restrict x < (γ : _)}
        add_mem' := by
          intro a b ha hb
          have ha' : Valued.v.restrict a < (γ : _) := ha
          have hb' : Valued.v.restrict b < (γ : _) := hb
          show Valued.v.restrict (a + b) < (γ : _)
          exact lt_of_le_of_lt (Valuation.map_add _ a b) (max_lt ha' hb')
        zero_mem' := by
          show Valued.v.restrict (0 : K) < (γ : _)
          simp [zero_lt_iff, Units.ne_zero γ]
        neg_mem' := by
          intro a ha
          have ha' : Valued.v.restrict a < (γ : _) := ha
          show Valued.v.restrict (-a) < (γ : _)
          simpa [Valuation.map_neg] using ha' }
    exact ⟨⟨H, H.isOpen_of_mem_nhds (g := 0) (Valued.mem_nhds_zero.2 ⟨γ, subset_rfl⟩)⟩, hγ⟩

/-- ★★★★**完備な付値環では「項が 0 に収束する」だけで総和可能**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★Tate 級数 `X(u,q) = Σ q^n u/(1−q^n u)² − 2s₁(q)` の収束はこれで判定できる。 -/
theorem summable_of_valued_tendsto_zero {α K Γ₀ : Type}
    [Ring K] [LinearOrderedCommGroupWithZero Γ₀] [Valued K Γ₀] [CompleteSpace K]
    (f : α → K) (hf : Tendsto f cofinite (𝓝 0)) : Summable f := by
  haveI := nonarchimedeanAddGroup_of_valued (K := K) (Γ₀ := Γ₀)
  exact (NonarchimedeanAddGroup.summable_iff_tendsto_cofinite_zero f).2 hf

/-! ## ★出典の紐付け(`.src`) -/

def nonarchimedeanAddGroup_of_valued.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数と局所高さ——Tate 級数の収束に要る非アルキメデス性)",
    sectionId := "genell-def-3-3" }

def summable_of_valued_tendsto_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数と局所高さ——Tate 級数の収束判定)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
