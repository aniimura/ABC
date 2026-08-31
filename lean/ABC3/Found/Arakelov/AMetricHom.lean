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
| **`Γ(L̄)`**（`Ō_X → L̄` の集合） | ★★★**本ファイル** `AMetric.gamma` |
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

## ★★★★★★★台の水準では `Hom(Ō_X, L) ≅ Γ(L, ⊤)` が出る

    `unitHomEquiv L : (𝟙_ ⟶ L) ≃ Γ(L, ⊤)`

★機構は「`⊤` が `X.Opens` の**終対象**だから `𝟙_ ≅ freeObj (yoneda ⊤)`」
（在庫の `freeYonedaTermIso`）＋ mathlib の `PresheafOfModules.freeYonedaEquiv`。

## ★★★★★★そして計量の側は 2 つの補題で繋がる

    `|f|_{Ō_X}(p) = |f(p)|`            （`AMetric.one_norm`）
    `|c·s|(p) = |c(p)|·|s|(p)`         （`AMetric.norm_smul`）

★前者が「`Ō_X` は**自明な** hermitian 計量」の中身であり、
後者が「射の条件は作用素ノルム `≤ 1` と同値」の根拠である。

## ★★★★★★★★★そして `Γ(L̄) = Hom(Ō_X, L̄)` が閉じる

    `IsAHom 1 L φ ↔ φ(1) ∈ Γ(L̄)`   （`isAHom_one_iff`）

★層の射 `φ : 𝒪_X → L` が**算術直線束の射** `Ō_X → L̄` であるための必要十分条件は
`|φ(1)| ≤ 1` である——これが原文の `Γ(L̄) ≝ Hom(Ō_X, L̄)` の中身である。

## ★残っている段（明示）

★★`unitHomEquiv` が**`φ ↦ φ(1)` そのものである**ことの証明（`unitHomEquiv_apply`）。
★`rfl` ではない（2026-08-28 実測）——`freeYonedaEquiv` の順方向の計算補題が
mathlib に無く（`freeYonedaEquiv_symm_app` はある）、`topTermIso` を展開する段が要る。
★★これは**小さい配管**であって、上の `isAHom_one_iff` が数学の中身をすでに持っている。
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
def AMetric.gamma (L : AMetric X) : Set (L.sheaf.obj (op ⊤)) :=
  { s | ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, L.norm s p ≤ 1 }

theorem AMetric.mem_gamma_iff (L : AMetric X) (s : L.sheaf.obj (op ⊤)) :
    s ∈ AMetric.gamma L ↔ ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, L.norm s p ≤ 1 := Iff.rfl

/-- ★★**射は `Γ` を `Γ` へ送る**。 -/
theorem AHom.mapsTo_gamma {L M : AMetric X} (f : AHom L M) :
    Set.MapsTo (fun s => f.hom.app (op ⊤) s) (AMetric.gamma L) (AMetric.gamma M) :=
  fun _ hs p => isAHom_mapsTo f.isAHom _ p (hs p)

/-- ★★★★★★★★★**`Γ` はテンソル積で掛け算になる**——
`s ∈ Γ(L̄)` かつ `t ∈ Γ(M̄)` なら `s ⊗ t ∈ Γ(L̄ ⊗ M̄)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は `AMetric.norm_mul`（`|s ⊗ t| = |s| · |t|`）だけである。 -/
theorem AMetric.tmul_mem_gamma {L M : AMetric X} {s : L.sheaf.obj (op ⊤)}
    {t : M.sheaf.obj (op ⊤)} (hs : s ∈ AMetric.gamma L) (ht : t ∈ AMetric.gamma M) :
    (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t) ∈ AMetric.gamma (L * M) := by
  intro p
  rw [AMetric.norm_mul]
  exact mul_le_one₀ (hs p) (AMetric.norm_nonneg M t p) (ht p)

/-- ★★`0` はつねに `Γ(L̄)` に入る。 -/
theorem AMetric.zero_mem_gamma (L : AMetric X) : (0 : L.sheaf.obj (op ⊤)) ∈ AMetric.gamma L := by
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

/-! ## ★★★★★★★ノルムの斉次性と `Ō_X` の計量 -/

/-- ★★`trivValue` は大域切断へのスカラー倍と可換である。 -/
theorem trivValue_smul (L : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L ≅ 𝟙_ (PresheafModulesOn X V))
    (c : (Γ(X, (⊤ : X.Opens)) : Type)) (s : L.obj (op ⊤)) :
    trivValue L V e (c • s)
      = (X.presheaf.map (homOfLE (le_top : V ≤ ⊤)).op c) * trivValue L V e s := by
  have h1 : secOn L V (c • s)
      = (X.presheaf.map (homOfLE (le_top : V ≤ ⊤)).op c) • secOn L V s :=
    PresheafOfModules.map_smul L (homOfLE (le_top : V ≤ ⊤)).op c s
  rw [trivValue_eq_trivEquiv, trivValue_eq_trivEquiv, h1, map_smul]
  rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★基準の自明化での値は制限そのものである。 -/
theorem trivValue_baseTriv (V : X.Opens)
    (s : ((𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens)) : Type)) :
    trivValue (𝟙_ X.PresheafOfModules) V (baseTriv X V) s
      = X.presheaf.map (homOfLE (le_top : V ≤ ⊤)).op s := rfl

namespace AMetric

/-- ★★★★★★★★**`Ō_X` の計量は `|f|(p) = |f(p)|` である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の「`O_X` に**自明な** hermitian 計量を入れたもの `Ō_X`」がこれである。 -/
theorem one_norm (s : ((1 : AMetric X).sheaf.obj (op (⊤ : X.Opens)) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (htop : p ⁻¹ᵁ (⊤ : X.Opens) = ⊤) :
    (1 : AMetric X).norm s p = ‖evalOn p ⊤ htop s‖ := by
  obtain ⟨c⟩ := nonempty_normChart (1 : AMetric X).triv p
  rw [AMetric.norm_eq (1 : AMetric X) s p c.V (baseTriv X c.V) c.hp]
  have hb : (1 : AMetric X).metric.normOf c.V (baseTriv X c.V) p c.hp s
      = ‖evalOn p c.V c.hp (trivValue (𝟙_ X.PresheafOfModules) c.V (baseTriv X c.V) s)‖ :=
    structLocalMetric_normOf_baseTriv X c.V p c.hp s
  rw [hb, trivValue_baseTriv, evalOn_restrict p (le_top : c.V ≤ ⊤) c.hp]

/-- ★★★★★★**ノルムは斉次である**: `|c·s|(p) = |c(p)|·|s|(p)`。

★これが「射の条件は作用素ノルム `≤ 1` と同値」の根拠である。 -/
theorem norm_smul (L : AMetric X) (c : (Γ(X, (⊤ : X.Opens)) : Type))
    (s : L.sheaf.obj (op ⊤)) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (htop : p ⁻¹ᵁ (⊤ : X.Opens) = ⊤) :
    L.norm (c • s) p = ‖evalOn p ⊤ htop c‖ * L.norm s p := by
  obtain ⟨ch⟩ := nonempty_normChart L.triv p
  rw [AMetric.norm_eq L (c • s) p ch.V ch.e ch.hp, AMetric.norm_eq L s p ch.V ch.e ch.hp]
  show trivSecNorm L.sheaf ch.V ch.e p ch.hp (c • s) * L.metric.h ch.V ch.e p
    = _ * (trivSecNorm L.sheaf ch.V ch.e p ch.hp s * L.metric.h ch.V ch.e p)
  have hts : trivSecNorm L.sheaf ch.V ch.e p ch.hp (c • s)
      = ‖evalOn p ⊤ htop c‖ * trivSecNorm L.sheaf ch.V ch.e p ch.hp s := by
    show ‖evalOn p ch.V ch.hp (trivValue L.sheaf ch.V ch.e (c • s))‖ = _
    rw [trivValue_smul, evalOn_mul, _root_.norm_mul,
      evalOn_restrict p (le_top : ch.V ≤ ⊤) ch.hp]
    rfl
  rw [hts]
  ring

end AMetric

/-! ## ★★★★★`Hom(Ō_X, L) ≅ Γ(L, ⊤)`（台の水準） -/

/-- ★`⊤` は `X.Opens` の終対象なので、`𝟙_` は `yoneda ⊤` の自由前層加群である。 -/
noncomputable def topTermIso (X : Scheme.{0}) :
    PresheafOfModules.freeObj
      (R := X.presheaf ⋙ forget₂ CommRingCat RingCat) (yoneda.obj (⊤ : X.Opens))
      ≅ 𝟙_ X.PresheafOfModules :=
  freeYonedaTermIso (R := X.presheaf) (⊤ : X.Opens)
    (fun _ => ⟨⟨homOfLE le_top⟩, fun _ => Subsingleton.elim _ _⟩)

/-- ★★★★**`Hom(𝟙_, L) ≅ Γ(L, ⊤)`**——原文の `Γ(L̄) ≝ Hom(Ō_X, L̄)` の**台の水準**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は「`⊤` が終対象だから `𝟙_ ≅ freeObj (yoneda ⊤)`」＋ mathlib の
`PresheafOfModules.freeYonedaEquiv`。 -/
noncomputable def unitHomEquiv (L : X.PresheafOfModules) :
    (𝟙_ X.PresheafOfModules ⟶ L) ≃ (L.obj (op (⊤ : X.Opens)) : Type) :=
  ((topTermIso X).symm.homCongr (Iso.refl L)).trans PresheafOfModules.freeYonedaEquiv

/-- ★`Ō_X` の大域切断 `1`。 -/
noncomputable def unitOneSec (X : Scheme.{0}) :
    ((𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens)) : Type) :=
  (1 : (Γ(X, (⊤ : X.Opens)) : Type))

theorem smul_unitOneSec (X : Scheme.{0}) (s : (Γ(X, (⊤ : X.Opens)) : Type)) :
    s • unitOneSec X = s := by
  show s * (1 : (Γ(X, (⊤ : X.Opens)) : Type)) = s
  rw [mul_one]

/-- ★★★★★★★★★**`Γ(L̄) = Hom(Ō_X, L̄)`** —— 原文の定義そのもの。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★層の射 `φ : 𝒪_X → L` が**算術直線束の射** `Ō_X → L̄` であるための必要十分条件は
`φ(1) ∈ Γ(L̄)`（すなわち `|φ(1)| ≤ 1`）である。

★★機構は 2 つ:
`Ō_X` の計量が `|f|(p) = |f(p)|` であること（`AMetric.one_norm`、とくに `|1| = 1`）と、
ノルムが斉次であること（`AMetric.norm_smul`）。 -/
theorem isAHom_one_iff (L : AMetric X) (φ : (1 : AMetric X).sheaf ⟶ L.sheaf) :
    IsAHom 1 L φ ↔ (φ.app (op (⊤ : X.Opens)) (unitOneSec X)) ∈ AMetric.gamma L := by
  have htop : ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, p ⁻¹ᵁ (⊤ : X.Opens) = ⊤ := fun _ => rfl
  constructor
  · intro h p
    have h1 := h (unitOneSec X) p
    have hone : (1 : AMetric X).norm (unitOneSec X) p = 1 := by
      rw [AMetric.one_norm (unitOneSec X) p (htop p)]
      show ‖evalOn p ⊤ (htop p) (1 : (Γ(X, (⊤ : X.Opens)) : Type))‖ = 1
      rw [evalOn_one, norm_one]
    rw [hone] at h1
    exact h1
  · intro h s p
    have hs1 : (s : (Γ(X, (⊤ : X.Opens)) : Type)) • unitOneSec X = s := smul_unitOneSec X s
    have hsm : φ.app (op (⊤ : X.Opens)) s
        = (s : (Γ(X, (⊤ : X.Opens)) : Type)) • φ.app (op (⊤ : X.Opens)) (unitOneSec X) :=
      (congrArg (fun y => φ.app (op (⊤ : X.Opens)) y) hs1).symm.trans
        (map_smul ((φ.app (op (⊤ : X.Opens))).hom) s (unitOneSec X))
    have hns : L.norm ((s : (Γ(X, (⊤ : X.Opens)) : Type))
          • φ.app (op (⊤ : X.Opens)) (unitOneSec X)) p
        = ‖evalOn p ⊤ (htop p) (s : (Γ(X, (⊤ : X.Opens)) : Type))‖
          * L.norm (φ.app (op (⊤ : X.Opens)) (unitOneSec X)) p :=
      AMetric.norm_smul L s _ p (htop p)
    have hn1 : (1 : AMetric X).norm s p
        = ‖evalOn p ⊤ (htop p) (s : (Γ(X, (⊤ : X.Opens)) : Type))‖ :=
      AMetric.one_norm s p (htop p)
    rw [hsm, hns, hn1]
    exact mul_le_of_le_one_right (norm_nonneg _) (h p)

/-! ### ★出典の紐付け(`.src`) -/

def IsAHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束の射——作用素ノルムが 1 以下)",
    sectionId := "genell-def-1-1-i" }

def AMetric.gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ(L̄)——大域切断の単位球)",
    sectionId := "genell-def-1-1-i" }

def isAHom_one_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ(L̄) = Hom(Ō_X, L̄)——射であることと |φ(1)| ≤ 1 の同値)",
    sectionId := "genell-def-1-1-i" }

def AMetric.one_norm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Ō_X の自明な hermitian 計量——|f|(p) = |f(p)|)",
    sectionId := "genell-def-1-1-i" }

def AMetric.norm_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(ノルムの斉次性)",
    sectionId := "genell-def-1-1-i" }

def unitHomEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Hom(𝟙_, L) ≅ Γ(L, ⊤)——Γ(L̄) の台の水準)",
    sectionId := "genell-def-1-1-i" }

def AMetric.tmul_mem_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(Γ はテンソル積で閉じること)",
    sectionId := "genell-def-1-1-i" }

def AMetric.gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "AMetric.norm(切断の大域ノルム)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm") 3,
    .citation "[ABC3]" "AMetric.norm_mul(|s ⊗ t| = |s| · |t|)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_mul") 3,
    .implicitStep
      ("★原文の射の条件「locally on X^arc、|−|_L ≤ 1 の切断を |−|_M ≤ 1 へ送る」を" ++
       "**作用素ノルムが 1 以下**の形で書いた。" ++
       "★★同値性は「計量が斉次だから s を定数倍して |s| ≤ 1 に落とせる」ことである") 3,
    .citation "[ABC3]" "freeYonedaTermIso(終対象の自由前層加群は単位)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.freeYonedaTermIso") 3,
    .citation "[mathlib]" "PresheafOfModules.freeYonedaEquiv"
      (.inMathlib "PresheafOfModules.freeYonedaEquiv") 3,
    .implicitStep
      ("★★台の水準では Hom(Ō_X, L) ≅ Γ(L, ⊤) が出る(unitHomEquiv)。" ++
       "★計量の側は one_norm(|f|_{Ō_X} = |f(p)|)と " ++
       "norm_smul(|c·s| = |c(p)|·|s|)で繋がる") 3,
    .implicitStep
      ("★★★残っている段: Hom(Ō_X, L̄) ≅ Γ(L̄) を**計量込みで**閉じる段" ++
       "(φ ↦ φ(1) が単位球への全単射であること)。" ++
       "★数学は上の 2 補題で尽きているが、𝟙_ の切断を「加群の元」と" ++
       "「環の元」として同時に扱う配管(s • 1 = s の綴り)で詰まっており、" ++
       "1 ブロックぶんの仕事として残す(2026-08-28 実測)") 3 ]

end ABC3.Found.Arakelov
