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

## ★★一般の場合に残るもの —— 構造を測った(2026-08-17 深夜)

★部品はすべて mathlib にある:

| 部品 | 内容 |
|---|---|
| `Scheme.Hom.resLE f W U e : U.toScheme ⟶ W.toScheme` | アフィン開集合への制限 |
| `resLE_comp_ι : f.resLE W U e ≫ W.ι = U.ι ≫ f` | 制限と包含の両立 |
| `resLE_app_top : (f.resLE W U e).app ⊤ = topIso.hom ≫ f.appLE W U e ≫ topIso.inv` | 切断の水準 |
| `comap_comp` | `comap I (f ≫ g) = (I.comap g).comap f` |
| `ideal_comap_of_isOpenImmersion` | 開埋め込みへの制限 |

★★筋道: `U.ι ≫ f = f.resLE W U e ≫ W.ι` の両辺に `comap` を当て、
`U.toScheme`, `W.toScheme` が**アフィンスキーム**であることから
本ファイルの `ideal_comap_eq_map_of_isAffine` を適用する。

## ★★★しかし本当の障害はそこではない

★★**任意のアフィン開集合 `U ⊆ X` に対し、`f(U)` を含むアフィン開集合 `W ⊆ Y`
が存在するとは限らない。**

★実際、本プロジェクトの用途(`x_F : Spec 𝓞_F ⟶ X`、`X` は非アフィン)では
像が 1 次元なので、1 つのアフィン開集合に収まらない。

★★★**したがって一般の場合には `U` を細かく取り直す(被覆する)段と、
そこから貼り合わせる段が要る。**

## ★★被覆と貼り合わせは**できる**(実測 2026-08-17 深夜)

★**貼り合わせ**: `Scheme.IdealSheafData.ext_of_iSup_eq_top` が mathlib にある——
アフィン開集合の族が `X` を覆い、各々で `ideal` が一致すれば 2 つの層は等しい。

★★**被覆**: 各点 `x : X` について、`f(x)` を含むアフィン開集合 `W ⊆ Y` を取り、
`f ⁻¹ᵁ W` の中に `x` を含むアフィン開集合 `U` を取ればよい。
**これは `lake env lean` で実際に構成して確認した**(`isBasis_affineOpens` 2 回)。

## ★★★残る本当の摩擦 —— `Γ` の同一視

`U.toScheme` の大域切断 `Γ(U.toScheme, ⊤)` と `Γ(X, U)` は
**同型だが等しくない**(`U.1.topIso`)。
★アフィンの場合の結果を `f.resLE W U h : U.toScheme ⟶ W.toScheme` に適用して
戻すには、この同型による共役を通す必要がある。

★部品は揃っている:
`resLE_app_top : (f.resLE W U e).app ⊤ = W.topIso.hom ≫ f.appLE W U e ≫ U.topIso.inv`
`ideal_comap_of_isOpenImmersion`(開埋め込みへの制限)。

★★★**数学の穴は無い。残るのは同型による共役の帳簿だけである。**

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
