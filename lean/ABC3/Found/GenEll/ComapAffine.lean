import ABC3.Found.GenEll.ComapLocal

/-!
# [GenEll] Proposition 1.4, (i) —— **アフィンでは `comap` は `Ideal.map` である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★テンソル積は要らなかった

`ComapLocal.lean` は「逆向きにはファイバー積の `Γ` の計算(`pullbackSpecIso`)が要る」
と書いた。★★★**それは誤りだった。** アフィンの場合は
**Galois 接続だけで**出る:

    `comap ⊣ map`(`Scheme.IdealSheafData.map_gc`)

★★機構:
`K ≝ ofIdealTop (Ideal.map φ I₀)` と置くと、
`I ≤ map K f` が `Ideal.le_comap_map`(`I₀ ≤ comap φ (map φ I₀)`)から出る。
Galois 接続を通すと `comap I f ≤ K`、すなわち逆向きの包含になる。

★★★**これで「正面から要ると思ったものが要らなかった」の 8 例目**になる。

## ★★得られるもの

アフィンでは `comap` が `Ideal.map` そのものになるので、
**`Ideal.map_mul` から `comap_mul` が出る**——
`Proposition 1.4, (i)`(高さの加法性)のアフィンの場合が閉じる。

## ★一般の場合に残るもの

`X`, `Y` が一般のスキームのときは、アフィン開集合への制限を経由する必要がある
(`comap_comp` + 開埋め込みへの制限 + `resLE`)。
★本ファイルは**アフィンの場合だけ**を取る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory CategoryTheory.Limits

variable {X Y : Scheme}

/-! ## ★★★逆向きの包含 —— Galois 接続で出る -/

set_option maxHeartbeats 1000000 in
/-- ★★★**アフィンでは `comap` は `Ideal.map` に含まれる**(逆向きの包含)。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**テンソル積を使わない。** `comap ⊣ map` の Galois 接続だけで出る。

★機構: `K ≝ ofIdealTop (Ideal.map φ I₀)` と置くと
`I ≤ map K f` が `Ideal.le_comap_map` から出て、
Galois 接続を通すと `comap I f ≤ K` になる。 -/
theorem ideal_comap_le_map_of_isAffine [IsAffine X] [IsAffine Y]
    (I : Y.IdealSheafData) (f : X ⟶ Y) :
    (I.comap f).ideal ⟨⊤, isAffineOpen_top X⟩
      ≤ (I.ideal ⟨⊤, isAffineOpen_top Y⟩).map (f.app ⊤).hom := by
  set J : Ideal Γ(X, ⊤) := (I.ideal ⟨⊤, isAffineOpen_top Y⟩).map (f.app ⊤).hom with hJ
  set K : X.IdealSheafData := Scheme.IdealSheafData.ofIdealTop J with hK
  have hKtop : K.ideal ⟨⊤, isAffineOpen_top X⟩ = J := by
    rw [hK]; simp
  have hle : I ≤ Scheme.IdealSheafData.map K f := by
    refine Scheme.IdealSheafData.le_of_isAffine ?_
    rw [Scheme.IdealSheafData.ideal_map_of_isAffineHom]
    have hKp : K.ideal ⟨f ⁻¹ᵁ (⊤ : Y.Opens), _⟩ = J := hKtop
    rw [hKp, hJ]
    exact fun x hx => Ideal.mem_comap.2 (Ideal.mem_map_of_mem _ hx)
  have h2 := (Scheme.IdealSheafData.le_map_iff_comap_le.1 hle) ⟨⊤, isAffineOpen_top X⟩
  rw [hKtop] at h2
  exact h2

/-! ## ★★★等式 -/

set_option maxHeartbeats 1000000 in
/-- ★★★**アフィンでは `comap` は `Ideal.map` である**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★易しい側は `ComapLocal.lean` の `map_appLE_le_ideal_comap`、
逆向きは上の `ideal_comap_le_map_of_isAffine`。

★★★**これで `comap` が完全に代数の言葉になった。** -/
theorem ideal_comap_eq_map_of_isAffine [IsAffine X] [IsAffine Y]
    (I : Y.IdealSheafData) (f : X ⟶ Y) :
    (I.comap f).ideal ⟨⊤, isAffineOpen_top X⟩
      = (I.ideal ⟨⊤, isAffineOpen_top Y⟩).map (f.app ⊤).hom := by
  refine le_antisymm (ideal_comap_le_map_of_isAffine I f) ?_
  have hle : (⊤ : X.Opens) ≤ f ⁻¹ᵁ (⊤ : Y.Opens) := le_top
  have h := map_appLE_le_ideal_comap I f ⟨⊤, isAffineOpen_top X⟩
    ⟨⊤, isAffineOpen_top Y⟩ hle
  have happ : f.appLE (⊤ : Y.Opens) (⊤ : X.Opens) hle = f.app ⊤ := by
    rw [Scheme.Hom.appLE]
    simp
  rwa [happ] at h

/-! ## ★★★`comap_mul`(アフィンの場合) -/

set_option maxHeartbeats 1000000 in
/-- ★★★**アフィンでは `comap` は積を保つ**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**これが `Proposition 1.4, (i)`(高さの加法性)の残り 1 本**である
——アフィンの場合について。

★機構は `ideal_comap_eq_map_of_isAffine` + `Ideal.map_mul` + `le_of_isAffine`。 -/
theorem comap_mul_of_isAffine [IsAffine X] [IsAffine Y]
    (I J : Y.IdealSheafData) (f : X ⟶ Y) :
    (I * J).comap f = I.comap f * J.comap f := by
  have key : ((I * J).comap f).ideal ⟨⊤, isAffineOpen_top X⟩
      = (I.comap f * J.comap f).ideal ⟨⊤, isAffineOpen_top X⟩ := by
    simp only [Scheme.IdealSheafData.ideal_mul, Pi.mul_apply,
      ideal_comap_eq_map_of_isAffine, Ideal.map_mul]
    rfl
  exact le_antisymm (Scheme.IdealSheafData.le_of_isAffine key.le)
    (Scheme.IdealSheafData.le_of_isAffine key.ge)

/-! ## ★出典の紐付け(`.src`)

★条つきである。一般のスキームではアフィン開集合への制限を経由する必要がある。 -/

def comap_mul_of_isAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(引き戻しが積を保つこと——アフィンの場合のみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
