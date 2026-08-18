import ABC3.Found.Arakelov.PicAwayLoc

/-!
# Arakelov (B1) 第 97 ブロック —— **`powers g` は `M_{t·g}` に可逆に作用する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★自由性を運ぶための可逆性

第 96 ブロックで `R_{t·g}` が `R_g` の局所化であることが出た。
★★加群の側で `M_g → M_{t·g}` を作るには、
**`powers g` が `M_{t·g}` に可逆に作用する**ことが要る。

★★★`g ∣ t·g` なので `g` は `R_{t·g}` で可逆であり、
その作用は**単元によるスカラー倍**だから全単射である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `bijective_smul_of_isUnit` | ★単元によるスカラー倍は全単射 |
| `isUnit_end_powers` | ★★★**`powers g` は `M_{t·g}` に可逆に作用する** |

## ★★★次

これで `IsLocalizedModule.lift` が使え、`M_g →ₗ M_{t·g}` が作れる。
★あとは `IsBaseChange.comp_iff`(§9-101)で局所化性を移し、
`Module.free_of_isLocalizedModule` で自由性が運べる。
-/

universe u

namespace ABC3.Found.Arakelov

/-- ★**単元によるスカラー倍は全単射である**。 -/
theorem bijective_smul_of_isUnit {A P : Type u} [CommRing A] [AddCommGroup P] [Module A P]
    {a : A} (h : IsUnit a) : Function.Bijective (fun x : P => a • x) := by
  obtain ⟨u, rfl⟩ := h
  have he : (fun x : P => (u : A) • x) = fun x => u • x := by
    funext x; rw [Units.smul_def]
  rw [he]
  exact (MulAction.toPerm u).bijective

variable (R : Type u) [CommRing R] (g t : R)
  (M : Type u) [AddCommGroup M] [Module R M]

/-- ★★★**`powers g` は `M_{t·g}` に可逆に作用する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`g ∣ t·g` だから `g` は `R_{t·g}` で可逆である。 -/
theorem isUnit_end_powers (s : Submonoid.powers g) :
    IsUnit (algebraMap R
      (Module.End R (LocalizedModule (Submonoid.powers (t * g)) M)) (s : R)) := by
  obtain ⟨n, hn⟩ := s.2
  rw [← hn, Module.End.isUnit_iff]
  have hu : IsUnit (algebraMap R (Localization (Submonoid.powers (t * g))) (g ^ n)) := by
    rw [map_pow]
    exact (IsLocalization.Away.isUnit_of_dvd (x := t * g) ⟨t, by ring⟩).pow n
  have hb := bijective_smul_of_isUnit
    (P := LocalizedModule (Submonoid.powers (t * g)) M) hu
  convert hb using 1
  funext x
  exact (IsScalarTower.algebraMap_smul (Localization (Submonoid.powers (t * g))) (g ^ n) x).symm

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_end_powers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——powers g は M_{t·g} に可逆に作用する)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
