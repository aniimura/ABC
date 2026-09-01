/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Meta.Claim

/-!
# 第 1241 ブロック —— **商の類の存在**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`quotLCyclicJ_spec` を使えるようにする

`quotLCyclicJ`（`EllModuliObjects`）は `∃ y, IsQuotClassJ x l y.1` のときだけ
意味を持つ（そうでなければ自分自身）。

☆位数 `l` の有理点 `Q` と、その Vélu の商が**楕円・半安定・乗法還元**であることから
その存在が出る。★個数の条件 `S.card + 1 = l` は第 1240。

## ★★逸脱の記録

☆商が**乗法還元を持つ**ことは仮説として受ける
（原文は「同種なので自動」と括弧で述べる。`VeluQuotOK` と同じ扱い）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep IsDedekindDomain NumberField ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★
**商の類は存在する**——★（第 1241）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆位数 `l` の有理点 `Q` と、その Vélu の商が楕円・半安定・乗法還元であることから
`IsQuotClassJ` の存在が出る。★個数の条件は第 1240。

★★逸脱: 商が乗法還元を持つことは仮説として受ける
（原文は「同種なので自動」と括弧で述べる）。 -/
theorem exists_isQuotClassJ (x : RealizedClass) {l : ℕ} (hl : l.Prime)
    (Q : x.rep.toSSCurve.W.toAffine.Point) (hQ : addOrderOf Q = l)
    (hell : (veluQuotientFull x.rep.toSSCurve.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
      SemistableAt p (veluQuotientFull x.rep.toSSCurve.W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))))
    (hmult : (quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))
      hell hss).HasMultRed) :
    ∃ y : RealizedClass, IsQuotClassJ x l y.1 := by
  refine ⟨⟨(quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).j,
    ⟨⟨quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss,
      hmult⟩, rfl⟩⟩, ?_⟩
  exact ⟨((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)), hell, hss,
    card_image_pointCoords_nsmul hl hQ, rfl⟩

/-! ## ★★★★★★★★★★★★商の類の `deg∞` -/

/-- ★★★★★★★★★★★★
**商の類の `deg∞` は `l` 倍**——★**無条件**（第 1242）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`degInfJ_eq`（類の値は曲線の `deg∞`）と
`degInfOf_eq_of_local`（局所の関係を足し上げる、在庫）を繋いだだけ。

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ` が要る形である
——局所の関係 `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` は `Lemma 3.2, (ii)`。 -/
theorem degInfJ_quot_eq (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ)
    (hloc : ∀ p : HeightOneSpectrum (𝓞 E.fld),
      minDeltaExp p (veluQuotientFull E.W S) = l * minDeltaExp p E.W) :
    degInfJ (quotSSCurve E S hell hss).j = (l : ℝ) * degInfJ E.j := by
  rw [degInfJ_eq, degInfJ_eq]
  exact ABC3.Found.GaloisRep.degInfOf_eq_of_local E.W (veluQuotientFull E.W S) l hloc

/-! ## ★★★★★★★★★★★★商の類の `ht^Falt` -/

/-- ★★★★★★★★★★★★
**商の類の `ht^Falt` の評価は曲線の評価そのもの**——★**無条件**（第 1243）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`faltingsHeightJ_eq`（類の値は曲線の `ht^Falt`、在庫）で置き換えるだけ。

★★★これが `EllModuliWitness` の `faltingsHeightJ_quotLCyclicJ` が要る形である
——同種写像の高さ評価 `hfalt` をそのまま類の言葉に移す。 -/
theorem faltingsHeightJ_quot_le (E : SSCurve) (S : Finset (E.fld × E.fld))
    (hell : (veluQuotientFull E.W S).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 E.fld), SemistableAt p (veluQuotientFull E.W S))
    (l : ℕ) (C₀ : ℝ)
    (hfalt : htFaltOf E.fld (veluQuotientFull E.W S)
      ≤ htFaltOf E.fld E.W + 2 * Real.log l + C₀) :
    faltingsHeightJ (quotSSCurve E S hell hss).j
      ≤ faltingsHeightJ E.j + 2 * Real.log l + C₀ := by
  rw [faltingsHeightJ_eq, faltingsHeightJ_eq]
  exact hfalt

/-! ## ★★★★★★★★★★★★★★★★`quotLCyclicJ` の水準へ -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `deg∞` は `l` 倍**——★（第 1244）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`quotLCyclicJ_spec`（在庫）で `S` を取り出し、第 1242 を当てる。
★局所の関係は `Lemma 3.2, (ii)` が与える（仮説として受ける）。

★★★これが `EllModuliWitness` の `degInfJ_quotLCyclicJ` そのものの形である。 -/
theorem degInfJ_quotLCyclicJ_of_local (x : RealizedClass) (l : ℕ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hloc : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
        minDeltaExp p (veluQuotientFull x.rep.toSSCurve.W S)
          = l * minDeltaExp p x.rep.toSSCurve.W) :
    degInfJ (quotLCyclicJ x l).cls = (l : ℝ) * degInfJ x.cls := by
  obtain ⟨S, hell, hss, hcard, hj⟩ := quotLCyclicJ_spec x l hex
  have hq := degInfJ_quot_eq x.rep.toSSCurve S hell hss l (hloc S hcard)
  rw [hj] at hq
  have hx : x.cls = x.rep.toSSCurve.j := (RealizedClass.rep_j x).symm
  rw [hx]
  exact hq

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `ht^Falt` の評価**——★（第 1245）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`quotLCyclicJ_spec`（在庫）で `S` を取り出し、第 1243 を当てる。
★同種写像の高さ評価は仮説として受ける。

★★★これが `EllModuliWitness` の `faltingsHeightJ_quotLCyclicJ` そのものの形である。 -/
theorem faltingsHeightJ_quotLCyclicJ_of_isog (x : RealizedClass) (l : ℕ) (C₀ : ℝ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hfalt : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      htFaltOf x.rep.toSSCurve.fld (veluQuotientFull x.rep.toSSCurve.W S)
        ≤ htFaltOf x.rep.toSSCurve.fld x.rep.toSSCurve.W + 2 * Real.log l + C₀) :
    faltingsHeightJ (quotLCyclicJ x l).cls
      ≤ faltingsHeightJ x.cls + 2 * Real.log l + C₀ := by
  obtain ⟨S, hell, hss, hcard, hj⟩ := quotLCyclicJ_spec x l hex
  have hq := faltingsHeightJ_quot_le x.rep.toSSCurve S hell hss l C₀ (hfalt S hcard)
  rw [hj] at hq
  have hx : x.cls = x.rep.toSSCurve.j := (RealizedClass.rep_j x).symm
  rw [hx]
  exact hq

/-! ## ★★★★★★★★★★★★局所の関係を全素点で -/

/-- ★★★★★★★★★★★★
**`Δ_min` の `l` 倍関係は全素点で成り立つ**——★**無条件**（第 1247）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆悪い素点（`v_p(j) < 0`）では `minDeltaExp_eq_mul_of_jExp_mul`（在庫）、
良い素点（`0 ≤ v_p(j)`）では両辺とも `0` である。

★★★これが第 1244（`quotLCyclicJ` の `deg∞`）が要る `hloc` の形であり、
残るのは `Lemma 3.2, (ii)` が与える `v_p(j(E′)) = l·v_p(j(E))` だけである。 -/
theorem minDeltaExp_eq_mul_of_jExp_all {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') {l : ℕ}
    (hbad : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E)
    (hgood : ∀ p : HeightOneSpectrum (𝓞 L), 0 ≤ jExp p E → 0 ≤ jExp p E') :
    ∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E := by
  intro p
  rcases lt_or_ge (jExp p E) 0 with hneg | hnn
  · exact minDeltaExp_eq_mul_of_jExp_mul p E E' (hss p) (hss' p) hneg (hbad p hneg)
  · have h1 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E (hss p)]
      exact max_eq_left (by omega)
    have h2 : minDeltaExp p E' = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E' (hss' p)]
      have hg := hgood p hnn
      exact max_eq_left (by omega)
    rw [h1, h2, mul_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `deg∞`（`jExp` の言葉で）**——★（第 1248）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1247（`Δ_min` の関係は `jExp` の関係から出る）を第 1244 に流し込んだ形。

★★★これで witness の `degInfJ_quotLCyclicJ` は
**`v_p(j(E′)) = l·v_p(j(E))`（`Lemma 3.2, (ii)`）ただ 1 つ**に帰着した。 -/
theorem degInfJ_quotLCyclicJ_of_jExp (x : RealizedClass) (l : ℕ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hrel : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      ∀ [_inst : (veluQuotientFull x.rep.toSSCurve.W S).IsElliptic],
      (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
        SemistableAt p (veluQuotientFull x.rep.toSSCurve.W S)) →
      (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
          jExp p x.rep.toSSCurve.W < 0 →
          jExp p (veluQuotientFull x.rep.toSSCurve.W S)
            = (l : ℤ) * jExp p x.rep.toSSCurve.W)
        ∧ (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
          0 ≤ jExp p x.rep.toSSCurve.W →
          0 ≤ jExp p (veluQuotientFull x.rep.toSSCurve.W S))) :
    degInfJ (quotLCyclicJ x l).cls = (l : ℝ) * degInfJ x.cls := by
  obtain ⟨S, hell, hss, hcard, hj⟩ := quotLCyclicJ_spec x l hex
  haveI := hell
  obtain ⟨hbad, hgood⟩ := hrel S hcard hss
  have hloc := minDeltaExp_eq_mul_of_jExp_all x.rep.toSSCurve.W
    (veluQuotientFull x.rep.toSSCurve.W S) x.rep.toSSCurve.ss hss hbad hgood
  have hq := degInfJ_quot_eq x.rep.toSSCurve S hell hss l hloc
  rw [hj] at hq
  have hx : x.cls = x.rep.toSSCurve.j := (RealizedClass.rep_j x).symm
  rw [hx]
  exact hq

/-- ★★★★★★★★★★★★
**悪い素点の関係を全素点へ束ねる**——★**無条件**（第 1258）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`hbad` は `minDeltaExp_eq_mul_at_bad_prime`
（`Skeleton/GenEll/TateIsogeny.lean`、証明済み）が与える形である。
★良い素点（`0 ≤ v_p(j)`）では両辺とも `0`。

★★★これが第 1244（`quotLCyclicJ` の `deg∞`）が要る `hloc` を
**局所の結果から直接**作る段である。 -/
theorem minDeltaExp_eq_mul_all {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') {l : ℕ}
    (hbad : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 →
      minDeltaExp p E' = l * minDeltaExp p E)
    (hgood : ∀ p : HeightOneSpectrum (𝓞 L), 0 ≤ jExp p E → 0 ≤ jExp p E') :
    ∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E := by
  intro p
  rcases lt_or_ge (jExp p E) 0 with hneg | hnn
  · exact hbad p hneg
  · have h1 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E (hss p)]
      exact max_eq_left (by omega)
    have h2 : minDeltaExp p E' = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E' (hss' p)]
      have hg := hgood p hnn
      exact max_eq_left (by omega)
    rw [h1, h2, mul_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `deg∞`（局所の `Δ_min` の言葉で）**——★（第 1259）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮説 `hbad` は `minDeltaExp_eq_mul_at_bad_prime`
（`Skeleton/GenEll/TateIsogeny.lean`、**証明済み**）が与える形そのものである。

★★★これで witness の `degInfJ_quotLCyclicJ` は
**在庫の局所結果 ＋ 良い素点で符号が保たれること**だけになった。 -/
theorem degInfJ_quotLCyclicJ_of_bad (x : RealizedClass) (l : ℕ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hrel : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      ∀ [_inst : (veluQuotientFull x.rep.toSSCurve.W S).IsElliptic],
      (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
        SemistableAt p (veluQuotientFull x.rep.toSSCurve.W S)) →
      (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
          jExp p x.rep.toSSCurve.W < 0 →
          minDeltaExp p (veluQuotientFull x.rep.toSSCurve.W S)
            = l * minDeltaExp p x.rep.toSSCurve.W)
        ∧ (∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
          0 ≤ jExp p x.rep.toSSCurve.W →
          0 ≤ jExp p (veluQuotientFull x.rep.toSSCurve.W S))) :
    degInfJ (quotLCyclicJ x l).cls = (l : ℝ) * degInfJ x.cls := by
  obtain ⟨S, hell, hss, hcard, hj⟩ := quotLCyclicJ_spec x l hex
  haveI := hell
  obtain ⟨hbad, hgood⟩ := hrel S hcard hss
  have hloc := minDeltaExp_eq_mul_all x.rep.toSSCurve.W
    (veluQuotientFull x.rep.toSSCurve.W S) x.rep.toSSCurve.ss hss hbad hgood
  have hq := degInfJ_quot_eq x.rep.toSSCurve S hell hss l hloc
  rw [hj] at hq
  have hx : x.cls = x.rep.toSSCurve.j := (RealizedClass.rep_j x).symm
  rw [hx]
  exact hq

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★
**`quotLCyclicJ` の `ht^Falt`（定数 `0` の形）**——★（第 1260）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮説は `htFalt_veluQuotientFull_le`
（`Found/GaloisRep/VeluNormalized.lean`、第 704、**証明済み**）が与える形そのもの
——定数は `0` でよい。

★★★これで witness の `faltingsHeightJ_quotLCyclicJ` は
**在庫の同種写像の高さ評価に幾何のデータを渡すだけ**になった。 -/
theorem faltingsHeightJ_quotLCyclicJ_zero (x : RealizedClass) (l : ℕ)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1)
    (hfalt : ∀ (S : Finset (x.rep.toSSCurve.fld × x.rep.toSSCurve.fld)),
      S.card + 1 = l →
      htFaltOf x.rep.toSSCurve.fld (veluQuotientFull x.rep.toSSCurve.W S)
        ≤ htFaltOf x.rep.toSSCurve.fld x.rep.toSSCurve.W + 2 * Real.log l) :
    faltingsHeightJ (quotLCyclicJ x l).cls
      ≤ faltingsHeightJ x.cls + 2 * Real.log l := by
  have h := faltingsHeightJ_quotLCyclicJ_of_isog x l 0 hex
    (fun S hS => by simpa using hfalt S hS)
  simpa using h

/-! ## ★★★★★★★★★★★★★★★★生成元を持ち歩く版 -/

open scoped Classical in
/-- ★★★★★★★★★★★★
**生成元 `Q` を持ち歩く `IsQuotClassJ`**（第 1262）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`IsQuotClassJ` は座標集合 `S` だけを持つので、
`htFalt_veluQuotientFull_le`（第 704）が要る**生成元 `Q`** が取り出せない。
★本定義は `Q` を持ち歩く形であり、`IsQuotClassJ` を含意する（下の定理）。

★★★これが `EllModuliWitness` の `faltingsHeight_quotLCyclic` を
在庫の同種写像の高さ評価に繋ぐための設計である。 -/
def IsQuotClassPointJ (x : RealizedClass) (l : ℕ) (y : ℂ) : Prop :=
  ∃ (Q : x.rep.toSSCurve.W.toAffine.Point) (_hQ : addOrderOf Q = l)
    (hell : (veluQuotientFull x.rep.toSSCurve.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic)
    (hss : ∀ p : HeightOneSpectrum (𝓞 x.rep.toSSCurve.fld),
      SemistableAt p (veluQuotientFull x.rep.toSSCurve.W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))),
    (quotSSCurve x.rep.toSSCurve
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).j = y

open scoped Classical in
/-- ★★★★★★★★**生成元つきなら座標集合つきでもある**——★**無条件**（第 1262）。

☆個数の条件 `S.card + 1 = l` は第 1240 が与える。 -/
theorem isQuotClassJ_of_point {x : RealizedClass} {l : ℕ} (hl : l.Prime) {y : ℂ}
    (h : IsQuotClassPointJ x l y) : IsQuotClassJ x l y := by
  obtain ⟨Q, hQ, hell, hss, hj⟩ := h
  exact ⟨((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)), hell, hss,
    card_image_pointCoords_nsmul hl hQ, hj⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def degInfJ_quotLCyclicJ_of_jExp.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の deg∞——jExp の言葉で。Lemma 3.2, (ii) を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def IsQuotClassPointJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(生成元 Q を持ち歩く IsQuotClassJ)",
    sectionId := "genell-lemma-3-5" }

def isQuotClassJ_of_point.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(生成元つきなら座標集合つきでもある。★無条件)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quotLCyclicJ_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の ht^Falt——定数 0 の形)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ_of_bad.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の deg∞——局所の Δ_min の言葉で)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_all.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(悪い素点の関係を全素点へ束ねる。★無条件)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_jExp_all.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Δ_min の l 倍関係は全素点で成り立つ。★無条件)",
    sectionId := "genell-lemma-3-2" }

def faltingsHeightJ_quotLCyclicJ_of_isog.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の ht^Falt の評価——同種写像の高さ評価を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(quotLCyclicJ の deg∞ は l 倍——局所の関係は Lemma 3.2, (ii) を仮説で受ける)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quot_le.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の ht^Falt の評価は曲線の評価そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quot_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の deg∞ は l 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類の存在——商の乗法還元は仮説として受ける)",
    sectionId := "genell-lemma-3-5" }

def exists_isQuotClassJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "card_image_pointCoords_nsmul(第 1240、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.card_image_pointCoords_nsmul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1241）**——`quotLCyclicJ_spec`（`EllModuliObjects`）を" ++
       "使えるようにする段である。☆★★逸脱: 商が**乗法還元を持つ**ことは" ++
       "仮説として受ける（原文は「同種なので自動」と括弧で述べる。" ++
       "`VeluQuotOK` と同じ扱い）。") 3 ]

end ABC3.Found.GenEll
