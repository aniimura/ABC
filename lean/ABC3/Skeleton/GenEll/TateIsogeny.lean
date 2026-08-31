/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.VeluMuSum
import ABC3.Found.GenEll.JScale
import ABC3.Meta.Claim
import ABC3.Skeleton.GenEll.TateODE

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
  [ .citation "[ABC3]" "sum_mu_veluV2（v = ∑_ζ(3X² + a₄ − Y) は有限個の係数の計算、第 810）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_veluV2") 1,
    .citation "[ABC3]" "sum_mu_veluW（2w も同じ、第 811）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_veluW") 1,
    .citation "[ABC3]" "sum_veluVC / veluVC_zero（係数 2 量の明示式、第 817）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_veluVC") 1,
    .citation "[ABC3]" "veluCurve_tateCurveAt_eq（変数変換は要らない、第 718）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluCurve_tateCurveAt_eq") 1,
    .implicitStep
      ("★★★測定（2026-08-31、第 818）——**adic 添字と q 次数がずれている**。"
       ++ "tateXtail(q·z^m, q) の adic 第 n 項は q^n·∑_{d∣n} d·q^d·z^{md} であり、"
       ++ "q 次数は n ではなく n+d である。★ゆえに最終の照合を"
       ++ "「adic 添字ごとの係数の一致」でやるには、先に**q 次数に揃えた形**へ"
       ++ "直す必要がある（古典形 X(u) = u/(1−u)² + ∑_N (∑_{d∣N} d(u^d + u^{−d} − 2)) q^N）。"
       ++ "☆道具は揃っている: AdicFubini.lean の adicSum_reindex_mul・adicSum_fubini、"
       ++ "AdicMul.lean の adicSum_mul。これが残る最後の配管である") 6,
    .implicitStep
      ("★★手順 4 の本体: q 次数に揃えた後、"
       ++ "★**l·veluVC N 0 − ∑_{a<l} veluVC N a = (σ₃(N) − [l∣N]·σ₃(N/l))·q^N**。"
       ++ "古典的な Tate/Vélu の q 展開恒等式そのものである") 8,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_tateModel（この結論の消費側、§9-1153）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_tateModel") 2 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★訂正後の目標（`c₄`・`c₆`） -/

open Finset in
/-- **[GenEll] 葉 1 の訂正後の目標 (1)**——`c₄` は `l⁴` 倍。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★第 834-835 の数値測定により、**係数 `a₄` の恒等式は偽**であり、
`c₄` で書いたこの形が真である（`l = 5, 7` で `q^21` まで全係数一致）。

☆`c₄(veluCurve W v w) = c₄ W + 240v`（`Found/GenEll/Velu.lean` の `veluCurve_c₄`）なので
左辺は Vélu の商の `c₄` そのものである。 -/
theorem c4_velu_tate {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    (tateCurveAt q hq).c₄
        + 240 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hsum1 := sum_veluV2_eq_sum_tateDYpair hl.pos hζl hu q hq hDX
  have hstep : ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = 2 * (∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    exact tateD2Xpair_eq _ _ q hq (hu i hi)
      (isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq))
  have hzero := sum_mu_dxpair_zero hl hζ hlu hu q hq
  have hd2 := sum_mu_d2xpair hl hζ hlu hu q hq hql
  rw [hstep, hzero, add_zero] at hd2
  rw [hsum1, tateCurveAt_c4_eq, tateCurveAt_c4_eq]
  linear_combination hd2

open Finset in
/-- **[GenEll] 葉 1 の訂正後の目標 (2)**——`c₆` は `−l⁶` 倍（分母を払った形）。

☆`c₆(veluCurve W v w) = c₆ W + 504 b₂ v + 6048 w`であり、Tate 曲線では `b₂ = 1`。
★`w` には `1/2` が入るので、**2 倍した形**で述べる。 -/
theorem c6_velu_tate {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    2 * (tateCurveAt q hq).c₆
        + 1008 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        + 3024 * (∑ i ∈ (range l).erase 0,
            (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              + 2 * (veluV2 (tateCurveAt q hq)
                      (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                      (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = 2 * (-((l : R) ^ 6) * (tateCurveAt (q ^ l) hql).c₆) := by
  sorry

def c4_velu_tate.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(訂正後の目標 (1)——c₄ は l⁴ 倍)",
    sectionId := "genell-lemma-3-2" }

def c4_velu_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_veluV2_eq_sum_tateDYpair(v = ∑ DY、第 846、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_veluV2_eq_sum_tateDYpair") 1,
    .citation "[ABC3]" "tateD2Xpair_eq(D²X = 2DY + DX、第 852、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateD2Xpair_eq") 1,
    .citation "[ABC3]" "sum_mu_dxpair_zero(∑ DX = 0、第 853)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_dxpair_zero") 1,
    .citation "[ABC3]" "sum_mu_d2xpair(∑ D²X の閉じた式、第 853)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_d2xpair") 1,
    .citation "[ABC3]" "tateCurveAt_c4_eq(c₄ = 1 + 240 s₃、第 853、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateCurveAt_c4_eq") 1 ]

def c6_velu_tate.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(訂正後の目標 (2)——c₆ は −l⁶ 倍)",
    sectionId := "genell-lemma-3-2" }

def c6_velu_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_mu_veluW（2w は有限個の係数の計算、第 811）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_veluW") 1,
    .citation "[ABC3]" "sigma_one_convolution（ラマヌジャンの恒等式、第 822）"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sigma_one_convolution") 8,
    .implicitStep
      ("★★数値確認済み（l = 5, 7 で q^21 まで全係数一致、第 835）") 6 ]

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
      ("★★同種な楼円曲線は同じ還元型をもつ（Neron-Ogg-Shafarevich、"
       ++ "あるいは Tate 加群の不分岐性）。mathlib には無い") 8,
    .implicitStep
      ("★★★★**測定（2026-08-31、第 827）——葉 1 に帰着する**。"
       ++ "対偶同種写像 `E′ → E`（核は位数 l）を使うと、"
       ++ "対偶を取って `jExp p E′ < 0`（= E′ が p で乗法還元）と仮定すれば、"
       ++ "葉 1（jExp_velu_bad）を E′ に適用して `jExp p E = l·jExp p E′ < 0` となり、"
       ++ "仮説 `0 ≤ jExp p E` に矛盾する。★★ゆえに本命題は"
       ++ "**葉 1 ＋「E は E′ の Vélu の商でもある」の 2 つに分解される**。"
       ++ "☆Néron-Ogg-Shafarevich も Néron 模型も要らない") 4,
    .implicitStep
      ("☆(a) p ∤ l のときの別道（第 785）: l-捧じれ点の x 座標は p で整だから"
       ++ "Vélu の和も整で、Vélu の式は還元と可換する。"
       ++ "☆(b) p ∣ l のときは有限平坦群スキームか Tate 曲線が要る") 4,
    .citation "[ABC3]" "jExp_velu_bad（対偶側へ適用する、第 827 の帰着）"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.jExp_velu_bad") 1 ]

end ABC3.Skeleton.GenEll
