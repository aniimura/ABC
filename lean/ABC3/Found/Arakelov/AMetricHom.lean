/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricNorm
import ABC3.Meta.Claim

/-!
# 算術直線束の**射**と `Γ(L̄)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★原文が `Definition 1.1, (i)` で挙げる 4 つ

| 原文 | 本プロジェクト |
|---|---|
| 算術直線束 `L̄ = (L, |−|_L)` | `AMetric`（`AMetricMonoid.lean`） |
| **射 `L̄ → M̄`**（`|−|_L ≤ 1` を `|−|_M ≤ 1` へ送る） | ★★**本ファイル** `IsAHom` / `AHom` |
| **`Γ(L̄)`**（`Ō_X → L̄` の集合） | ★★★**本ファイル** `arithGamma` |
| テンソル積と `APic(X)` | `AMetricPic.lean` |

原文の射の条件は

> such that locally on `X^arc`, sections of `L` with `|−|_L ≤ 1`
> map to sections of `M` with `|−|_M ≤ 1`

である。★これは「作用素ノルムが `1` 以下」——すなわち
`|φ(s)|_M(p) ≤ |s|_L(p)` ——と同値であり、本ファイルはその形で書く。
★★同値性は「計量が斉次だから、`s` を定数倍して `|s| ≤ 1` に落とせる」こと。

## ★★★★★★`Γ(L̄)` は**単位球**である

原文は `Γ(L̄) ≝ Hom(Ō_X, L̄)` と定める。
★`Ō_X` の切断 `f` は `|f|(p) = |f(p)|` を持ち、`f ↦ f·s` の像は `|f·s| = |f|·|s|`
なので、条件は `|s|_L ≤ 1` に他ならない。★★したがって

    `Γ(L̄) = { s ∈ Γ(L) | ∀ p, |s|_L(p) ≤ 1 }`

## ★残っている段（明示）

★`Hom(Ō_X, L̄) ≅ Γ(L̄)` を**射の側から**証明する段（`Γ(L) ≅ Hom(𝒪_X, L)` の
`⊤` の水準の形）。★★在庫は制限した圏（`PresheafModulesOn X V`）の水準にある
（`unitHomOfSection`、`PicTrivialNoSheaf.lean`）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★★★★★射 -/

/-- ★★★★★★**算術直線束の射であること**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の「`|−|_L ≤ 1` の切断を `|−|_M ≤ 1` の切断へ送る」を
**作用素ノルムが `1` 以下**の形で書いたもの。 -/
def IsAHom (L M : AMetric X) (φ : L.sheaf ⟶ M.sheaf) : Prop :=
  ∀ (s : L.sheaf.obj (op ⊤)) (p : Spec (CommRingCat.of ℂ) ⟶ X),
    M.norm (φ.app (op ⊤) s) p ≤ L.norm s p

/-- ★★**射は単位球を単位球へ送る**——原文の条件そのもの。 -/
theorem isAHom_mapsTo {L M : AMetric X} {φ : L.sheaf ⟶ M.sheaf} (h : IsAHom L M φ)
    (s : L.sheaf.obj (op ⊤)) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (hs : L.norm s p ≤ 1) : M.norm (φ.app (op ⊤) s) p ≤ 1 :=
  le_trans (h s p) hs

/-- ★★★**算術直線束の射**。 -/
structure AHom (L M : AMetric X) where
  /-- 台となる層の射。 -/
  hom : L.sheaf ⟶ M.sheaf
  /-- ★ノルムを増やさないこと。 -/
  isAHom : IsAHom L M hom

namespace AHom

/-- ★恒等射。 -/
def id (L : AMetric X) : AHom L L :=
  ⟨𝟙 L.sheaf, fun s p => le_of_eq (congrArg (fun x => L.norm x p) rfl)⟩

/-- ★★合成。 -/
def comp {L M N : AMetric X} (f : AHom L M) (g : AHom M N) : AHom L N where
  hom := f.hom ≫ g.hom
  isAHom := fun s p => by
    have h1 := g.isAHom (f.hom.app (op ⊤) s) p
    have h2 := f.isAHom s p
    have hcomp : (f.hom ≫ g.hom).app (op ⊤) s
        = g.hom.app (op ⊤) (f.hom.app (op ⊤) s) := rfl
    rw [hcomp]
    exact le_trans h1 h2

end AHom

/-! ## ★★★★★★`Γ(L̄)` -/

/-- ★★★★★★★★**`Γ(L̄)` —— 大域切断の単位球**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文は `Γ(L̄) ≝ Hom(Ō_X, L̄)` と定めるが、
`Ō_X` の切断 `f` の像は `f·s` で `|f·s| = |f|·|s|` なので、
条件は `|s|_L ≤ 1` に他ならない。 -/
def arithGamma (L : AMetric X) : Set (L.sheaf.obj (op ⊤)) :=
  { s | ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, L.norm s p ≤ 1 }

theorem mem_arithGamma_iff (L : AMetric X) (s : L.sheaf.obj (op ⊤)) :
    s ∈ arithGamma L ↔ ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, L.norm s p ≤ 1 := Iff.rfl

/-- ★★**射は `Γ` を `Γ` へ送る**。 -/
theorem AHom.mapsTo_arithGamma {L M : AMetric X} (f : AHom L M) :
    Set.MapsTo (fun s => f.hom.app (op ⊤) s) (arithGamma L) (arithGamma M) :=
  fun _ hs p => isAHom_mapsTo f.isAHom _ p (hs p)

/-- ★★★★★★★★★**`Γ` はテンソル積で掛け算になる**——
`s ∈ Γ(L̄)` かつ `t ∈ Γ(M̄)` なら `s ⊗ t ∈ Γ(L̄ ⊗ M̄)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は `AMetric.norm_mul`（`|s ⊗ t| = |s| · |t|`）だけである。 -/
theorem tmul_mem_arithGamma {L M : AMetric X} {s : L.sheaf.obj (op ⊤)}
    {t : M.sheaf.obj (op ⊤)} (hs : s ∈ arithGamma L) (ht : t ∈ arithGamma M) :
    (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t) ∈ arithGamma (L * M) := by
  intro p
  rw [AMetric.norm_mul]
  exact mul_le_one₀ (hs p) (AMetric.norm_nonneg M t p) (ht p)

/-- ★★`0` はつねに `Γ(L̄)` に入る。 -/
theorem zero_mem_arithGamma (L : AMetric X) : (0 : L.sheaf.obj (op ⊤)) ∈ arithGamma L := by
  intro p
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  rw [AMetric.norm_eq L 0 p c.V c.e c.hp]
  show trivSecNorm L.sheaf c.V c.e p c.hp 0 * L.metric.h c.V c.e p ≤ 1
  have hev : evalOn p c.V c.hp 0 = 0 := by
    show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      (((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq c.hp.symm)).op).hom
        ((p.app c.V).hom 0)) = 0
    rw [map_zero, map_zero, map_zero]
  have h0 : trivSecNorm L.sheaf c.V c.e p c.hp 0 = 0 := by
    show ‖evalOn p c.V c.hp (trivValue L.sheaf c.V c.e 0)‖ = 0
    rw [trivValue_zero, hev, norm_zero]
  rw [h0, zero_mul]
  exact zero_le_one

/-! ### ★出典の紐付け(`.src`) -/

def IsAHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束の射——作用素ノルムが 1 以下)",
    sectionId := "genell-def-1-1-i" }

def arithGamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ(L̄)——大域切断の単位球)",
    sectionId := "genell-def-1-1-i" }

def tmul_mem_arithGamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ はテンソル積で閉じること)",
    sectionId := "genell-def-1-1-i" }

def arithGamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "AMetric.norm(切断の大域ノルム)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm") 3,
    .citation "[ABC3]" "AMetric.norm_mul(|s ⊗ t| = |s| · |t|)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_mul") 3,
    .implicitStep
      ("★原文の射の条件「locally on X^arc、|−|_L ≤ 1 の切断を |−|_M ≤ 1 へ送る」を" ++
       "**作用素ノルムが 1 以下**の形で書いた。" ++
       "★★同値性は「計量が斉次だから s を定数倍して |s| ≤ 1 に落とせる」ことである") 3,
    .implicitStep
      ("★★★残っている段: Hom(Ō_X, L̄) ≅ Γ(L̄) を射の側から証明する段" ++
       "(Γ(L) ≅ Hom(𝒪_X, L) の ⊤ の水準の形)。" ++
       "★在庫は制限した圏(PresheafModulesOn X V)の水準にある(unitHomOfSection)") 3 ]

end ABC3.Found.Arakelov
