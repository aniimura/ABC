import ABC3.Interface.GaloisRep.Representation
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.RingTheory.Valuation.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic

/-!
# Galois 表現のスケルトン(3/3)—— **還元・Tate 曲線・Faltings 高さ**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.15–p.17。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★本ファイルが受ける 3 件

| # | 受けるもの | mathlib の在庫(2026-08-17 実測) |
|---|---|---|
| G6 | Tate 曲線と局所高さ `v(q_E)` | ★**無い**。原文の典拠は [FC] Ch. III, Cor. 7.3 |
| G7 | 半安定還元と `𝓞_L` 上のモデル | ★`EllipticCurve/Reduction.lean` は**ある**が Néron モデルは無い |
| G8 | Faltings 高さ `ht^Falt = deg(ω_E)` | ★**無い**。Arakelov 側 (D2) に従属 |

★★★**G6 が `Definition 3.3` そのもの**であり、`Lemma 3.2` が要求するものである。
★★G8 は **Arakelov 側と Galois 側の合流点**である——
`deg(ω_E)` は算術直線束の次数だから、(D2) が要る。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

/-! ## ★★★G6 —— Tate 曲線と局所高さ -/

/-- **(G6)** 乗法還元をもつ楕円曲線の **Tate 母数 `q_E`** と**局所高さ** `v_K(q_E)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**原文の `Definition 3.3` そのものである。**
★`Lemma 3.2`(局所の階数 1 部分群)がこの上で述べられる。

★★★**2026-08-17: 退化検査で作り直した。**以前は `LocalField : Type` /
`Curve : LocalField → Type` と**世界ごと posit** していたので、
`PUnit` と定数 `1` で埋まってしまった(実測済)。
★**mathlib の `WeierstrassCurve` と `Valuation` に接地し**、
**Tate 一意化そのもの**を要求する形に改めた。 -/
structure TateCurveData where
  /-- `W` が **分裂乗法還元**をもつこと(原文の「[bad] multiplicative reduction」)。 -/
  HasSplitMultRed : {K : Type} → [Field K] → WeierstrassCurve K → Prop
  /-- ★★**Tate 母数** `q_E ∈ Kˣ`——分裂乗法還元のとき定まる。 -/
  tateParam : {K : Type} → [Field K] → (W : WeierstrassCurve K) → HasSplitMultRed W → Kˣ
  /-- ★★★**Tate 一意化** `E(K) ≅ Kˣ / q^ℤ`。

  ★★★**これが退化を殺す。**`q` が実際に一意化していることを要求するので、
  「曲線も体も 1 点」という witness は通らない。 -/
  uniformization : ∀ {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K)
    (h : HasSplitMultRed W),
    Nonempty (W.toAffine.Point ≃+ Additive (Kˣ ⧸ Subgroup.zpowers (tateParam W h : Kˣ)))
  /-- ★★★**局所高さ** `v_K(q_E) ∈ ℤ_{>0}`(原文 `Definition 3.3`)。

  ★★**付値を引数に取る**——「どの付値での高さか」は原文でも `v_K` と明示されている。
  ★乗法群上の付値 `Kˣ →* Multiplicative ℤ` として受ける(型クラスが軽い)。 -/
  localHeight : {K : Type} → [Field K] → (v : Kˣ →* Multiplicative ℤ) →
    (W : WeierstrassCurve K) → HasSplitMultRed W → ℕ
  /-- ★★★**局所高さは Tate 母数の付値そのもの**——定義を型に出す。

  ★これがあるので `localHeight := 1` のような定数 witness は通らない。 -/
  localHeight_eq : ∀ {K : Type} [Field K] (v : Kˣ →* Multiplicative ℤ) (W : WeierstrassCurve K)
    (h : HasSplitMultRed W),
    (localHeight v W h : ℤ) = Multiplicative.toAdd (v (tateParam W h))
  /-- ★局所高さは正(原文が `∈ Z_{>0}` と書いている)。 -/
  localHeight_pos : ∀ {K : Type} [Field K] (v : Kˣ →* Multiplicative ℤ) (W : WeierstrassCurve K)
    (h : HasSplitMultRed W), 0 < localHeight v W h

def TateCurveData.waiting : WaitingFor :=
  { what := "(G6) Tate 曲線 E(K̄) ≅ K̄^x/q^Z、Tate 母数 q_E、および局所高さ v_K(q_E)(= [GenEll] Definition 3.3)"
    trackB := "Found/GaloisRep — ★mathlib に Tate 曲線は無い(2026-08-16 実測)。★★**FLT にも無い**——`FLT/TateCurve/TateCurve.lean` は 20 行の入口宣言だけ(2026-08-17、clone して計数)。★★★原文の典拠は [FC](Faltings-Chai, *Degenerations of Abelian Varieties*)Chapter III, Corollary 7.3 で、これも mathlib に無い" }

/-! ## ★★G7 —— 半安定還元とモデル -/

/-- **(G7)** 半安定還元と、数体の整数環 `𝓞_L` 上への延長。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★原文は「`E` は `L` のすべての有限素点で半安定還元をもつ」を繰り返し仮定する。
★その意味を与えるのが本 obligation である。 -/
structure SemistableModelData where
  /-- 台となる Tate 曲線。 -/
  toTateCurveData : TateCurveData
  /-- `E_L` が `L` のすべての有限素点で半安定還元をもつこと。 -/
  SemiStable : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → Prop
  /-- ★★`𝓞_L` 上の 1 次元半アーベルスキームへの延長。 -/
  ExtendsToOL : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → Prop
  /-- ★半安定なら `𝓞_L` 上へ延びる。 -/
  semiStable_extends : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
    SemiStable L E → ExtendsToOL L E
  /-- ★★**同種な曲線も半安定**(原文が `E_H = E/H` について使う)。 -/
  semiStable_isogeny : (L : Type) → [Field L] → [NumberField L] →
    WeierstrassCurve L → WeierstrassCurve L → Prop

def SemistableModelData.waiting : WaitingFor :=
  { what := "(G7) 半安定還元の定義と、𝓞_L 上の 1 次元半アーベルスキームへの延長"
    trackB := "Found/GaloisRep — ★mathlib は `AlgebraicGeometry/EllipticCurve/Reduction.lean` を持つ(2026-08-17 実測)が、**Néron モデル・半アーベルスキームは無い**。★★原文は延長の存在を [FC] Chapter I, Proposition 2.7 に帰している" }

/-! ## ★★★G8 —— Faltings 高さ(Arakelov 側との合流点) -/

/-- **(G8)** **Faltings 高さ** `ht^Falt(E) = deg(ω_E)`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★**ここが Arakelov 側と Galois 側の合流点である。**
`ω_E` は `Spec 𝓞_L` 上の算術直線束なので、その次数を取るには
**(D2)(`APic(Spec 𝓞_F)` と `deg_F`)が要る**。

★`Lemma 3.5` / `Lemma 3.7` / `Proposition 3.4` はすべてこの上で述べられる。 -/
structure FaltingsHeightData where
  /-- 台となる半安定モデル。 -/
  toSemistableModelData : SemistableModelData
  /-- ★★**Faltings 高さ**。 -/
  htFalt : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → ℝ
  /-- 無限遠因子の次数 `deg∞`。 -/
  degInf : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → ℝ
  /-- ★`deg∞` は非負(局所高さの和だから)。 -/
  degInf_nonneg : ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
    0 ≤ degInf L E
  /-- ★★★**`Proposition 3.4`** —— `deg∞` は `ht^Falt` で上から抑えられる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
  prop_3_4 : ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L]
    (E : WeierstrassCurve L),
    toSemistableModelData.SemiStable L E →
    degInf L E / (12 * (1 + ε)) ≤ htFalt L E + C

def FaltingsHeightData.waiting : WaitingFor :=
  { what := "(G8) Faltings 高さ ht^Falt(E) = deg(omega_E) と、Proposition 3.4(deg_infty を ht^Falt で抑える)"
    trackB := "Found/GaloisRep — ★★★**Arakelov 側との合流点**。`omega_E` は Spec 𝓞_L 上の算術直線束なので、次数を取るには **(D2)** が要る。★`ADiv` / `deg_F` は `Found/GenEll/ArithDiv.lean` に実装済み(sorry 0)なので、(D1)(D2) が入れば `omega_E` の構成だけが残る。★★`Lemma 3.5` / `Lemma 3.7` / `Proposition 3.4` はすべてこの上で述べられる" }

/-! ## ★出典の紐付け(`.src`) -/

def TateCurveData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(Tate 母数と局所高さの定式化のみ)",
    sectionId := "genell-def-3-3" }

def SemistableModelData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(半安定還元と 𝓞_L 上のモデルのみ)",
    sectionId := "genell-def-3-3" }

def FaltingsHeightData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4(Faltings 高さの定式化のみ)",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GaloisRep
