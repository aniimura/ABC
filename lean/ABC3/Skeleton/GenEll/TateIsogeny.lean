/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Meta.Claim

/-!
# `Lemma 3.2, (ii)` の曲線の水準 —— **`E_q/μ_l` は母数 `q^l` の Tate 曲線**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か

`Found/GaloisRep/DegInfLocal.lean` の `minDeltaExp_eq_mul_of_tateModel`（`§9-1153`、第 716）は

> `E` の Tate モデルの母数が `q`、`E′` の母数が `q^l` なら
> `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`

を**無条件で**証明している。★したがって `Lemma 3.5`（および `Lemma 3.7` の第 3 の主張）に
残っている入力は、**本ファイルが固定する 1 つの命題だけ**である:

> **`E′ = E/μ_l` は母数 `q^l` の Tate モデルを持つ。**

## ★★★なぜ `Skeleton` に置くか

`Lemma 3.2` は本プロジェクトで**完了扱い**の項目だが、`Found/GenEll/Lemma32.lean` の
`lemma_3_2` が取っているのは**群論の核**

    `(Lˣ/q^ℤ) ⧸ (x ↦ x^l の核の像) ≃* Lˣ/(q^l)^ℤ`

であって、**曲線の水準の主張**（`E_q/μ_l` の Weierstrass モデルが `tateCurveAt (q^l)` と
変数変換で一致すること）ではない。★★この差が `Lemma 3.5` の最後の穴である。

## ★★★★★証明の筋（`q`-展開の恒等式）

`tateCurveAt q` は `⟨1, 0, 0, tateA4(q), tateA6(q)⟩`（`tateA4 = −5·σ₃`）である。
`veluCurve` は `a₁, a₂, a₃` を変えないので、**変数変換は要らず係数の一致を見ればよい**:

    `tateA4(q^l) = tateA4(q) − 5·v`,
    `tateA6(q^l) = tateA6(q) − b₂·v − 7·w`   （`b₂ = a₁² + 4a₂ = 1`）

ここで `v = Σ_{ζ ∈ μ_l∖{1}} g^x(P_ζ)`、`w = Σ (u_{P_ζ}/2 + g^x(P_ζ)·x_{P_ζ})` であり、
`P_ζ` の座標は `tateXpair ζ (qζ⁻¹) q`・`tateYpair ζ (qζ⁻¹) q`（`Found/GaloisRep/TateOrigin.lean`）。

★★`ζ` について足すと `ζ` の指数が `l` の倍数の項だけが残る——これが
`σₖ(q) → σₖ(q^l)` を生む機構である。

☆本プロジェクトの Tate 冪級数の在庫（`Found/GaloisRep/Tate*.lean` が 97 ファイル）で
計算できる。**新しい数学は要らない。**

## ★次に何をすれば終わりか

1. `μ_l` の点の座標を `q`-冪級数として書く（`tateXpair`・`tateYpair` の特殊化）
2. `Σ_{ζ ∈ μ_l}` を取ると `ζ` の指数が `l` の倍数の項だけ残ることを示す
3. 上の 2 式（`a₄`・`a₆`）を確かめる
4. `minDeltaExp_eq_mul_of_tateModel`（第 716）へ渡す
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-- **[GenEll] `Lemma 3.2, (ii)` の曲線の水準**——`E_q/μ_l` は母数 `q^l` の Tate 曲線。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★★★★★☆**これが `Lemma 3.5` に残る唯一の穴である。**
`minDeltaExp_eq_mul_of_tateModel`（`§9-1153`、第 716）がこの結論を消費して
`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` を出す。

☆仮説 `hquot` は「`E′` は `E` を位数 `l` の巡回部分群で割ったもので、
その部分群は各悪い素点で `μ_l` に対応する」という `Lemma 3.5` の設定
（原文の *global rank one subgroup*）を型にしたものである。 -/
theorem tateModel_of_quot_mu {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W W' : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    [W'.IsElliptic] [W'.IsMinimal R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (l : ℕ) (hl : 0 < l)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (D : VariableChange R) (hD : D • integralModel R W = tateCurveAt q hq)
    (hquot : True) :
    ∃ D' : VariableChange R,
      D' • integralModel R W' = tateCurveAt (q ^ l) hql := by
  sorry

def tateModel_of_quot_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(曲線の水準——E_q/μ_l は母数 q^l の Tate 曲線)",
    sectionId := "genell-lemma-3-2" }

def tateModel_of_quot_mu.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★`q`-展開の恒等式: tateA4(q^l) = tateA4(q) − 5v、" ++
       "tateA6(q^l) = tateA6(q) − v − 7w（b₂ = 1）。" ++
       "v・w は μ_l の点にわたる Vélu の和である") 8,
    .citation "[ABC3]" "sum_mu_tateXpair_eq（∑_ζ X(ζ,q) は定数項を除いて ζ-free、第 795）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_tateXpair_eq") 1,
    .citation "[ABC3]" "sum_mu_tateYpair_eq（∑_ζ Y(ζ,q) も同じ、第 796）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_tateYpair_eq") 1,
    .citation "[ABC3]" "twelve_mul_sum_mu_ringInverse（12·∑ ζ/(1−ζ)² = −(l²−1)、第 794）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.twelve_mul_sum_mu_ringInverse") 1,
    .citation "[ABC3]" "sum_mu_frac_cube（∑ ζ²/(1−ζ)³ = (l²−1)/24、第 797）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_frac_cube") 1,
    .citation "[ABC3]" "sum_mu_adicSum_mul（積を μ_l 上で足す道具、第 798）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_adicSum_mul") 1,
    .implicitStep
      ("★★★進捗（2026-08-31、第 786-798）: 手順 1・2は**済んだ**。" ++
       "μ_l 上の指標和（MuCharSum.lean、24 本）と " ++
       "I 進和と有限和の可換性（AdicFinsetSum.lean）により、" ++
       "∑_ζ X(ζ,q)・∑_ζ Y(ζ,q) はどちらも `[l ∣ d]·l − 1` を係数に持つ " ++
       "adicSum になった。定数項も両方計算済み。" ++
       "☆残るのは手順 3——v = ∑(3x² + a₄ − y)・w の組み立てと a₄・a₆ の照合。" ++
       "★そのためには ζ-次数を mod l で揃えた「μ-等級付き I 進級数」の枠が要る") 8,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_tateModel（この結論の消費側、§9-1153）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_tateModel") 2 ]

/-! ## ★★★★★★★★★★★★★★★★`Lemma 3.5` が直接消費する形（`j` の付値） -/

/-- **[GenEll] `Lemma 3.5` の残る入力 (1)**——悪い還元の素点での `v_p(j′) = l·v_p(j)`。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★`Found/GaloisRep/Lemma35Concrete.lean` の `lemma_3_5_velu_j`（`§9-1177`、第 751）が
**そのまま消費する形**である。

☆機構: 悪い還元の素点では `E` は Tate 曲線 `E_q` であり `v_p(j) = −v_p(q)`。
`H` が `μ_l` に対応するとき `E′ = E_q/μ_l = E_{q^l}` なので
`v_p(j′) = −l·v_p(q) = l·v_p(j)`。★これが `tateModel_of_quot_mu` の内容である。

☆仮説 `hmu` は「`H` は各悪い素点で `μ_l` に対応する」という `Lemma 3.5` の設定
（原文の *global rank one subgroup*）を型にしたものである。 -/
theorem jExp_velu_bad {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : Nat.Prime l) (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hmu : True)
    (p : HeightOneSpectrum (𝓞 L)) (hss : SemistableAt p E) (hbad : jExp p E < 0) :
    jExp p E' = (l : ℤ) * jExp p E := by
  sorry

/-- **[GenEll] `Lemma 3.5` の残る入力 (2)**——良い還元は同種写像で保たれる。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`v_p(j) ≥ 0` かつ半安定なら `E` は `p` で良い還元をもつ。★★同種な曲線は
同じ還元型をもつ（Néron–Ogg–Shafarevich、あるいは Tate 加群の不分岐性）ので
`E′` も良い還元をもち、`v_p(j′) ≥ 0` である。

☆こちらは `v_p(j′) = l·v_p(j)` **より弱く**てよい——半安定なら
`v_p(Δ_min) = max(0, −v_p(j))` なので、良い還元の素点では両辺とも `0` になる。 -/
theorem jExp_velu_good {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : Nat.Prime l) (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (p : HeightOneSpectrum (𝓞 L)) (hss : SemistableAt p E) (hss' : SemistableAt p E')
    (hgood : 0 ≤ jExp p E) :
    0 ≤ jExp p E' := by
  sorry

def jExp_velu_bad.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力 (1)——悪い素点での v_p(j′) = l·v_p(j))",
    sectionId := "genell-lemma-3-5" }

def jExp_velu_bad.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateModel_of_quot_mu(E_q/μ_l は母数 q^l の Tate 曲線)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateModel_of_quot_mu") 15,
    .implicitStep
      ("★Tate 曲線では v_p(j) = −v_p(q)(TateJ.lean)。母数が q^l なら " ++
       "v_p(j′) = −l·v_p(q) = l·v_p(j)") 4 ]

def jExp_velu_good.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力 (2)——良い還元は同種写像で保たれる)",
    sectionId := "genell-lemma-3-5" }

def jExp_velu_good.needs : List ProofObligation :=
  [ .implicitStep
      ("★★同種な楕円曲線は同じ還元型をもつ(Néron-Ogg-Shafarevich、" ++
       "あるいは Tate 加群の不分岐性)。mathlib には無い") 8 ,
    .implicitStep
      ("\u2605\u2605\u2605\u5206\u89e3(2026-08-31\u3001\u7b2c 785): \u672c\u547d\u984c\u306f\u7d20\u70b9\u3092 2 \u3064\u306b\u5206\u3051\u308b\u3068" ++
       "\u7247\u65b9\u306f Neron-Ogg-Shafarevich \u3092\u8981\u3089\u306a\u3044\u3002" ++
       "\u2606(a) p \u2224 l \u306e\u3068\u304d: l-\u6367\u3058\u308c\u70b9\u306e x \u5ea7\u6a19\u306f\u5206\u5272\u591a\u9805\u5f0f " ++
       "\u03c8_l(\u5148\u982d\u4fc2\u6570 l)\u306e\u6839\u3060\u304b\u3089 p \u3067\u6574\u3001\u3086\u3048\u306b Velu \u306e\u548c v, w \u3082\u6574\u3067\u3001" ++
       "Velu \u306e\u5f0f\u306f\u9084\u5143\u3068\u53ef\u63db\u3059\u308b\u3002\u2605E \u304c\u826f\u3044\u9084\u5143\u3092\u3082\u3064\u306a\u3089\u5546\u3082\u975e\u7279\u7570\u306a\u306e\u3067 " ++
       "v_p(\u0394') = 0\u3001\u3059\u306a\u308f\u3061 v_p(j') \u2265 0\u3002" ++
       "\u2606(b) p \u2223 l \u306e\u3068\u304d: \u6709\u9650\u5e73\u5766\u7fa4\u30b9\u30ad\u30fc\u30e0\u304b Tate \u66f2\u7dda\u304c\u8981\u308b\u3002" ++
       "\u2605\u2605(a) \u306f\u9084\u5143\u306e\u6a5f\u69cb\u3060\u3051\u3067\u6e08\u3080\u306e\u3067\u3001\u5148\u306b\u9589\u3058\u3089\u308c\u308b\u5074\u3067\u3042\u308b") 4,
    .implicitStep
      ("\u2606\u5225\u9053(\u6e2c\u5b9a\u6e08\u307f\u3001\u7b2c 785): \u30e2\u30b8\u30e5\u30e9\u30fc\u591a\u9805\u5f0f \u03a6_l(X,Y) \u2208 \u2124[X,Y] \u304c " ++
       "Y \u306b\u3064\u3044\u3066\u30e2\u30cb\u30c3\u30af\u3067 \u03a6_l(j, j') = 0 \u306a\u306e\u3067\u3001" ++
       "v_p(j) \u2265 0 \u306a\u3089 j' \u3082 p \u3067\u6574\u2014\u2014\u3053\u308c\u306f**\u7d14\u7c8b\u306b\u4ee3\u6570\u7684**\u3067 Neron \u6a21\u578b\u3092\u8981\u3089\u306a\u3044\u3002" ++
       "\u2606\u305f\u3060\u3057 \u03a6_l \u306e\u6574\u6570\u6027\u30fb\u30e2\u30cb\u30c3\u30af\u6027\u3092\u51fa\u3059\u306b\u306f X_0(l) \u306e\u30e2\u30b8\u30e5\u30e9\u30fc\u95a2\u6570\u8ad6\u304c\u8981\u308b\u306e\u3067\u3001" ++
       "\u73fe\u6642\u70b9\u3067\u306f (a) \u306e\u9053\u306e\u65b9\u304c\u77ed\u3044") 2 ]

end ABC3.Skeleton.GenEll
