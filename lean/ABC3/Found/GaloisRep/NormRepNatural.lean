/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateClassMap

/-!
# 正規化代表の**自然性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★`galois-equivariant-tate-uniformization` の段 2 の核

`ResearchPaper/mathlib-gap.json` の `galois-equivariant-tate-uniformization` は
5 段に割ってある（2026-08-27）:

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | 固定した完備体での**選ばれた**同型 `Kˣ/q^ℤ ≃+ E_q(K)` | ★済（`tatePhiAddEquivAll`） |
| 2 | 体の写像に対する同変性 | ★本ファイルがその**核**を取る |
| 3 | 有限拡大に対する自然性 | 未 |
| 4 | `K̄` への帰納極限 | 未 |
| 5 | `l`-捩れの完全列 | ★消費側は済（`Lemma32StableLine`） |

## ★★★★★★なぜ核なのか

`tatePhi` は `tateAOf` / `tateWOf` を通して定義されており、それらは
**`normRep`（正規化代表）で特徴づけられている**（`tateAOf_spec`）。
★したがって同変性は **`normRep` の自然性**に帰着する。

★★`normRep` は `.choose` で定義されているが、
**一意性補題 `eq_normRep` がある**——
「類 `c` の代表で付値が `[0, v q)` に入るものは `normRep c` に一致する」。
★★★**一意なものは自動的に自然である**——これが本ファイルの 3 行の証明である。

## ★残っている段（明示）

★本ファイルは `normRep` の自然性だけを取る。
★★`tatePhi` 全体の同変性にはさらに `tatePtPair` の自然性が要り、
そこは**級数と `σ` を交換させる段**（`σ` の連続性＝付値の拡大の一意性）を含む。
-/

namespace ABC3.Found.GaloisRep

/-! ## ★一般の付値保存準同型に沿った自然性 -/

/-- ★★★★★★**正規化代表は付値を保つ群準同型で自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

> `σ (normRep c) = normRep (σ_* c)`

★★機構は**一意性 1 本**である（`eq_normRep`）——
`σ (normRep c)` は `σ_* c` の代表であり、付値は `σ` で保たれるので
条件 `[0, v q)` を満たす。したがって `normRep (σ_* c)` に一致する。

★★★**一意なものは自動的に自然である。** -/
theorem normRep_map {K L : Type} [Field K] [Field L]
    (v : Kˣ →* Multiplicative ℤ) (w : Lˣ →* Multiplicative ℤ)
    (q : Kˣ) (hq : 0 < vAdd v q) (σ : Kˣ →* Lˣ)
    (hQ : 0 < vAdd w (σ q))
    (hσv : ∀ u : Kˣ, vAdd w (σ u) = vAdd v u)
    (hle : Subgroup.zpowers q ≤ Subgroup.comap σ (Subgroup.zpowers (σ q)))
    (c : Kˣ ⧸ Subgroup.zpowers q) :
    σ (normRep v q hq c)
      = normRep w (σ q) hQ (QuotientGroup.map _ _ σ hle c) := by
  refine eq_normRep w (σ q) hQ _ (σ (normRep v q hq c)) ?_ ?_ ?_
  · show QuotientGroup.mk (σ (normRep v q hq c)) = _
    rw [← QuotientGroup.map_mk' (Subgroup.zpowers q) (Subgroup.zpowers (σ q)) σ hle]
    exact congrArg _ (normRep_mk v q hq c)
  · rw [hσv]; exact normRep_nonneg v q hq c
  · rw [hσv, hσv]; exact normRep_lt v q hq c

/-! ## ★★Galois 作用の形（`σ q = q`） -/

/-- ★`σ q = q` なら `zpowers q` は `comap σ (zpowers q)` に入る。 -/
theorem zpowers_le_comap_self {K : Type} [Field K] (q : Kˣ) (σ : Kˣ →* Kˣ) (hσq : σ q = q) :
    Subgroup.zpowers q ≤ Subgroup.comap σ (Subgroup.zpowers q) := by
  intro x hx
  obtain ⟨n, rfl⟩ := hx
  refine Subgroup.mem_comap.2 ?_
  rw [map_zpow, hσq]
  exact Subgroup.zpow_mem _ (Subgroup.mem_zpowers q) n

/-- ★★★★★★★**Galois 作用の形の自然性**（`σ q = q`、付値を保つ）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★これが `G_K` が `Kˣ/q^ℤ` に作用するときの同変性である。
★`q ∈ K` は `G_K` で固定され、付値は一意なので `σ` は等距離写像になる
——その 2 つが仮説 `hσq`・`hσv` である。 -/
theorem normRep_map_self {K : Type} [Field K]
    (v : Kˣ →* Multiplicative ℤ) (q : Kˣ) (hq : 0 < vAdd v q) (σ : Kˣ →* Kˣ)
    (hσq : σ q = q) (hσv : ∀ u : Kˣ, vAdd v (σ u) = vAdd v u)
    (c : Kˣ ⧸ Subgroup.zpowers q) :
    σ (normRep v q hq c)
      = normRep v q hq (QuotientGroup.map _ _ σ (zpowers_le_comap_self q σ hσq) c) := by
  refine eq_normRep v q hq _ (σ (normRep v q hq c)) ?_ ?_ ?_
  · show QuotientGroup.mk (σ (normRep v q hq c)) = _
    rw [← QuotientGroup.map_mk' (Subgroup.zpowers q) (Subgroup.zpowers q) σ
      (zpowers_le_comap_self q σ hσq)]
    exact congrArg _ (normRep_mk v q hq c)
  · rw [hσv]; exact normRep_nonneg v q hq c
  · rw [hσv]; exact normRep_lt v q hq c

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**同変性の段の核**であって、
`tatePhi` 全体の同変性ではない。 -/

def normRep_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——正規化代表の自然性)",
    sectionId := "genell-def-3-3" }

def normRep_map_self.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Galois 作用に対する正規化代表の同変性)",
    sectionId := "genell-def-3-3" }

def normRep_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_normRep(正規化代表の一意性)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.eq_normRep") 15,
    .implicitStep
      ("★★★★**一意なものは自動的に自然である**——" ++
       "normRep は .choose で定義されているが、eq_normRep が一意性を与えるので" ++
       "σ (normRep c) が条件を満たすことを見れば済む") 15,
    .implicitStep
      ("★★残る段: tatePhi 全体の同変性には tatePtPair の自然性が要り、" ++
       "そこは**級数と σ を交換させる段**(σ の連続性 = 付値の拡大の一意性)を含む。" ++
       "★mathlib-gap の galois-equivariant-tate-uniformization の段 2 の残りである") 15,
    .implicitStep
      ("★逸脱: σ を「付値を保つ群準同型 Kˣ →* Lˣ」として受けている。" ++
       "★★体の自己同型から来る場合、hσv(付値保存)は" ++
       "付値の拡大の一意性から出るが、本ファイルはそれを仮説に置いている") 15 ]

end ABC3.Found.GaloisRep
