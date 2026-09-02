/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GenEll.QuotClassExists
import ABC3.Found.GenEll.Lemma37Hdag
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.DetCycloChar
import ABC3.Interface.GenEll.EllModuli
import ABC3.Skeleton.GenEll.Section3
import ABC3.Skeleton.GenEll.TateLocalModelK
import ABC3.Skeleton.GenEll.VeluSemistable
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Skeleton.GenEll.GaloisImage
import ABC3.Skeleton.GenEll.Section4
import ABC3.Meta.Claim

/-!
# `EllModuliData` の witness（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.23。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`Found/GenEll/EllModuliObjects.lean`・`EllModuliGalois.lean` で作った対象を
**実際に `EllModuliData` に組み上げる**。★これで 50 個近い欄の型が一度に検査される。

## ★★★残る `sorry` は 4 本の葉に対応する 6 欄だけである

| 欄 | 依存する葉 |
|---|---|
| `degInf_quotLCyclic` | 葉 1・2（`Lemma 3.5`） |
| `faltingsHeight_quotLCyclic` | 葉 1・2 |
| `galoisFinite_lcyclicExc` | 葉 1・2（`Lemma 3.5` の (†)）＋ `northcott`（済） |
| `imageContainsSL2_of_torsionExt` | 葉 3（局所理論の行列表示） |
| `imageSurjective_of_containsSL2` | 葉 4（円分指標の全射性） |

☆`mem_lcyclicExc` は `lcyclicExc` を**その定義そのもの**に取るので自明に埋まる
——有限性の側（`galoisFinite_lcyclicExc`）が本体である。

## ★★★★★★★★★★★★★★★★残りの地図（2026-09-02、第 1239-1249）

☆**4 本の `sorry` はすべて `Found` 側の形に降りた**。

| # | `sorry` | `Found` 側の形 | 残る前提 |
|---|---|---|---|
| 1 | `galoisFiniteJ_lcyclicExcJ` | `exists_galoisFinite_lcyclic`（第 1239） | `VeluQuotOK` |
| 2 | `degInfJ_quotLCyclicJ` | `degInfJ_quotLCyclicJ_of_jExp`（第 1248） | `Lemma 3.2, (ii)` |
| 3 | `faltingsHeightJ_quotLCyclicJ` | `faltingsHeightJ_quotLCyclicJ_of_isog`（第 1245） | 同種写像の高さ評価 |
| 4 | `imageContainsSL2J_torsionExt` | `imageContainsSL2J_of_galTate'`（第 1249） | `σ` の幂単性・非自明性 |

★★**測定**——(4) の `hclosed`（像の閉性）は**不要だった**。
`galRep` の連続性は `galRep_continuous'`（第 772、葉 5）で既に閉じている。

☆(2) の道筋: `h4`・`h6`（Vélu の `v, w` と Eisenstein の関係）
→ `j_velu_tate_eq_map`（在庫）→ `jExp_eq_mul_of_j_tate_pow`（第 997）
→ `minDeltaExp_eq_mul_of_jExp_all`（第 1247）→ `degInfJ_quotLCyclicJ_of_jExp`（第 1248）。
★★☆**訂正（2026-09-02、第 1256-1257）**——`h4`・`h6` は
`Skeleton/GenEll/TateIsogeny.lean` の `c4_velu_tate`・`c6_velu_tate` として
**既に証明済み**である（sorry 0）。
同じく (3) の「同種写像の高さ評価」も
`htFalt_veluQuotientFull_le`（第 704、証明済み）である。

★★★**したがって 4 本に残るのは配管だけであり、
未証明の外部引用は 1 本もない**。

## ★★★★★★★★★★★★★★★★★★★★**2026-09-02（第 1341-1342）——#1 が閉じた**

☆第 1341——界面の `galoisFinite_lcyclicExc` の `C` の量化を
`∀ C` から `∃ C₀, ∀ C ≥ C₀` に直した（`∀ C` は**作れない主張**だった）。
★第 1342——`Skeleton/GenEll/VeluSemistable.lean` に
**残った 1 本を節点として立てた**（`veluQuotOK_semistable`）。
☆`veluQuotOK_all` が `VeluQuotOK` を与えるので、
`galoisFiniteJ_lcyclicExcJ`（#1）は**本ファイルの `sorry` ではなくなった**。

★★★**本ファイルの `sorry` は 2 本**（#2 の「商の類の存在」と #4）であり、
別に `Skeleton/GenEll/VeluSemistable.lean` に名前のついた葉が 1 本ある。
☆これは**数が減ったのではなく、何が残っているかが 1 行になった**と読むこと。

## ★★★★★★★★★★★★★★★★★★★★**2026-09-02（第 1340）——#2 の残りは「存在」だけになった**

☆界面の `degInf_quotLCyclic` を**不等式に弱めた**（消費側はその向きしか使わない）。
★`degInfJ_quotLCyclicJ_ge_of_exists`（下）で**商の類が存在すれば無条件**に出る。

### ★★★これで消えたもの——良い素点の半安定性

等式のままなら良い素点で `minDeltaExp p E′ = 0`
（＝**同種で良還元が保たれる**＝Néron–Ogg–Shafarevich）が要った。
☆不等式なら `0 ≤ minDeltaExp p E′` で済むので、**その要求は消えた**。

### ★残っているものの測定（2026-09-02）

| # | 残る一点 | 性質 |
|---|---|---|
| 1 | `VeluQuotOK` の**良い素点の半安定性** | ★既知数学（Néron–Ogg–Shafarevich） |
| 2 | `HasLCyclicJ` から `IsQuotClassJ` の**存在** | ☆配管＋#1 と同じ半安定性 |
| 4 | 分裂乗法還元の準備と `ζ_l ∈ L_v` | ☆局所の配管 |

★★**#1 と #2 は同じ根**である——`quotSSCurve` が `SSCurve` を作る以上、
商が半安定であることは**構造上避けられない**（`degInfJ`・`faltingsHeightJ` が
`SSCurve` の代表で定義されているため）。
☆悪い素点側は第 1327、楕円性は第 1336 で閉じているので、
残るのは**「良い素点で Vélu の商が良還元を持つ」1 本**である。
★その根は剰余体の上の Vélu の定理（`Δ(Ẽ/H̃) ≠ 0`）であり、
本プロジェクトが `ℂ` で使った一意化の道（第 1330-1335）は標数 `p` では使えない。

## ★★★★★★★★★★★★★★★★★★★★**2026-09-02（第 1338-1339）——#3 が閉じた**

☆`faltingsHeightJ_quotLCyclicJ`（#3）は **`sorry` ではなくなった**。道は 2 段:

1. 第 1337-1338——`htFalt_veluQuotientFull_le_uncond` を**無条件**にした
   （一意化の族は選択で取れ、`hfin` は第 1083・1149 が無条件に与える）。
2. 第 1339——`IsQuotClassJ` を**生成元 `Q` を持ち歩く形**に締め直した
   （かつては任意の座標集合 `S` で、同種写像の評価を当てられなかった）。

★★★**残る `sorry` は 3 本**（#1 `galoisFiniteJ_lcyclicExcJ`・#2 `degInfJ_quotLCyclicJ`・#4 `imageContainsSL2J_torsionExt`）である。
☆#1・#2 が共通に欲しがるのは「Vélu の商の半安定性」（良い素点側）であり、
悪い素点側は第 1327、楕円性は第 1336 で閉じている。

## ★★★★★★★★★★★★★★★★★★★★この 4 本が 3 項目を同時に閉じる（第 1261）

☆**測定（2026-09-02）**——`Skeleton/GenEll/Section4.lean` は **sorry 0**
（`cor_4_3`・`cor_4_4` は `EllModuliData` から導出済み、第 367）であり、
`Skeleton/GenEll/GaloisImage.lean` の `theorem_3_8` も sorry 0 である。

★★★したがって**本ファイルの 4 本を閉じれば、
`Theorem 3.8`（§3 9/9）と `Corollary 4.3`・`4.4`（§4 5/5）が同時に落ちる**
——指標は 20/24 → **23/24** になる。
☆残るのは §2 の `Theorem 2.1`（本論文外の 2 結果）だけである。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll ABC3.Found.GaloisRep

open scoped Classical in
/-- ★★★★**`lcyclicExc` 欄**——`mem_lcyclicExc` が要求する条件そのものを集めた集合。

★有限性（`galoisFinite_lcyclicExc`）が本体であり、それは `Lemma 3.5` の (†) と
`northcott`（`§9-1179`、第 753）から出る。 -/
def lcyclicExcJ (C eps : ℝ) (KV : Set ℂ) : Set ℂ :=
  {x : ℂ | ∃ (E : RealizedClass) (l : ℕ), x = E.cls ∧ Nat.Prime l ∧
    E.rep.toSSCurve.SemiStable ∧ HasLCyclicJ E.rep.toSSCurve l ∧
    E.rep.toSSCurve.PrimeToLocalHeights l ∧
    (((100 * (E.degOfDefinition : ℝ)
          * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ eps) ≤ (l : ℝ))
        ∧ E.rep.toSSCurve.HasMultRed)
      ∨ E.cls ∈ KV)}

/-! ## ★★★★★★★★★★★★★★★★`lcyclicExc` の有限性

★★☆**測定の履歴**——第 1263 で「条件 (a) に `HasMultRed` が無いと
Galois-有限にならない」と測った。☆**その直しは適用済み**である
（上の `lcyclicExcJ` の条件 (a) に `HasMultRed` が入っている）。

★★★**2026-09-02（第 1341）**——さらに界面の `C` の量化を
`∀ C` から `∃ C₀, ∀ C ≥ C₀` に直した。☆`C` をいくらでも小さく取れると
条件 (a) が常に満たされてしまい、**作れない主張**だったからである。
-/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`lcyclicExcJ` は `VeluQuotOK` さえあれば `Galois`-有限**——★（第 1341）。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆`exists_galoisFinite_lcyclic`（第 1239）が与える `Exc` に**部分集合として入る**
（`galoisFiniteJ_subset`、第 1238）。

★★★これで #1 の残りは **`VeluQuotOK` ただ 1 つ**になった
——第 1336 によりそれは「Vélu の商の半安定性」だけである。 -/
theorem galoisFiniteJ_lcyclicExcJ_of_veluQuotOK
    (hquot : ∀ (E : SSCurve) (l : ℕ), VeluQuotOK E l)
    (eps : ℝ) (heps : 0 < eps) (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) :
    ∃ C₀ : ℝ, ∀ C : ℝ, C₀ ≤ C → GaloisFiniteJ (lcyclicExcJ C eps KV) := by
  obtain ⟨C₀, -, Exc, hExc, hmem⟩ := exists_galoisFinite_lcyclic KV hKV eps heps hquot
  refine ⟨C₀, fun C hC => galoisFiniteJ_subset ?_ hExc⟩
  intro x hx
  obtain ⟨E, l, rfl, hl, hss, hcyc, hpr, hcond⟩ := hx
  have hj : E.rep.j = E.cls := RealizedClass.rep_j E
  have hd0 : (0 : ℝ) ≤ (E.rep.deg : ℝ) := Nat.cast_nonneg _
  have hP : (0 : ℝ) ≤ (E.rep.deg : ℝ) ^ eps := Real.rpow_nonneg hd0 eps
  rw [← hj]
  refine hmem E.rep l hl ?_ hcyc
  rw [hj]
  rcases hcond with ⟨hle, hmult⟩ | hKVmem
  · have hdd : E.degOfDefinition = E.rep.deg := rfl
    rw [hdd] at hle
    refine Or.inl ⟨le_trans ?_ hle, hmult⟩
    have h1 : C₀ * (E.rep.deg : ℝ) ^ eps ≤ C * (E.rep.deg : ℝ) ^ eps :=
      mul_le_mul_of_nonneg_right hC hP
    nlinarith [hd0, h1]
  · exact Or.inr ⟨hKVmem, hpr⟩

/-- ★★★★★★★★★★★★★★★★
**`lcyclicExc` は `Galois`-finite**——★**これはもう `sorry` ではない**（2026-09-02、第 1342）。

☆`VeluQuotOK` は `Skeleton/GenEll/VeluSemistable.lean` の `veluQuotOK_all` から出る。
★その節点（`veluQuotOK_semistable`）は**良い素点の半安定性 1 本**であり、
本ファイルの `sorry` ではなくなった。

★★☆**測定（2026-09-02、第 1341）**——`galoisFiniteJ_lcyclicExcJ_of_veluQuotOK`（上）で
**`VeluQuotOK` を仮定すれば閉じる**ことを示した。
☆したがって本 `sorry` の中身は `VeluQuotOK`（第 1336 により「Vélu の商の半安定性」）
**ただ 1 つ**である。 -/
theorem galoisFiniteJ_lcyclicExcJ (eps : ℝ) (heps : 0 < eps) (KV : Set ℂ)
    (hKV : CompactlyBoundedJ KV) :
    ∃ C₀ : ℝ, ∀ C : ℝ, C₀ ≤ C → GaloisFiniteJ (lcyclicExcJ C eps KV) :=
  galoisFiniteJ_lcyclicExcJ_of_veluQuotOK veluQuotOK_all eps heps KV hKV

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**商の類が存在すれば `l·deg∞(E) ≤ deg∞(E′)`**——★**無条件**（第 1340）。

☆悪い素点では Tate の関係（`isMuAtBadPrimes_of_veluQuotient_nodeg`、第 1141）、
良い素点では `minDeltaExp p E = 0` なので自明である。

★★★**これが「不等式に弱める」ことの利得である**——等式なら
良い素点で `minDeltaExp p E′ = 0`（同種で良還元が保たれること）が要った。 -/
theorem degInfJ_quotLCyclicJ_ge_of_exists (x : RealizedClass) {l : ℕ} (hl : Nat.Prime l)
    (hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1) :
    (l : ℝ) * degInfJ x.cls ≤ degInfJ (quotLCyclicJ x l).cls := by
  obtain ⟨E, hEj, hEpr, Q, hQ, hell, hss, hj⟩ := quotLCyclicJ_spec x l hex
  haveI := hell
  have hcop := E.not_dvd_jExp_of_primeToLocalHeights hl hEpr
  have hmu := isMuAtBadPrimes_of_veluQuotient_nodeg E.W
    (veluQuotientFull E.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    hl Q hQ rfl E.ss hss hcop
  have hge := ABC3.Found.GaloisRep.degInfOf_ge_of_local E.W
    (veluQuotientFull E.W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) l
    (fun p => ABC3.Found.GaloisRep.minDeltaExp_le_of_bad_delta p _ _
      (E.ss p) l (fun hb => hmu p hb))
  have hq : (l : ℝ) * degInfJ E.j
      ≤ degInfJ (quotSSCurve E
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) hell hss).j := by
    rw [degInfJ_eq, degInfJ_eq]
    exact hge
  rw [hj, hEj] at hq
  exact hq

/-- ★★★`l·deg∞(E) ≤ deg∞(E′)`——★**残るのは「商の類の存在」だけ**（第 1340）。

★★☆**測定（2026-09-02、第 1340）**——不等式に弱めたことで
**良い素点の半安定性（同種で良還元が保たれること）は不要になった**。
☆残る `sorry` は「`HasLCyclicJ` から `IsQuotClassJ` の存在を出す」一点であり、
それは `Gal`-安定な直線の生成元を `L(H)` で有理化して降ろす段（第 1194-1199）である。 -/
theorem degInfJ_quotLCyclicJ (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hcyc : HasLCyclicJ x.rep.toSSCurve l)
    (hpr : x.rep.toSSCurve.PrimeToLocalHeights l) :
    (l : ℝ) * degInfJ x.cls ≤ degInfJ (quotLCyclicJ x l).cls := by
  by_cases hex : ∃ y : RealizedClass, IsQuotClassJ x l y.1
  · exact degInfJ_quotLCyclicJ_ge_of_exists x hl hex
  · sorry

/-- ★★★★★★★★★★★★★★★★★★★★
**`ht^Falt(E′) ≤ ht^Falt(E) + 2log(l)`——★これはもう `sorry` ではない**（2026-09-02、第 1339）。

☆第 1338（無条件の同種写像の高さ評価）と、`IsQuotClassJ` を生成元 `Q` を持ち歩く形に締め直したことで、
定数 `C₀ = 0` で閉じた。 -/
theorem faltingsHeightJ_quotLCyclicJ :
    ∃ C₀ : ℝ, ∀ (x : RealizedClass) (l : ℕ), Nat.Prime l → HasLCyclicJ x.rep.toSSCurve l →
      faltingsHeightJ (quotLCyclicJ x l).cls
        ≤ faltingsHeightJ x.cls + 2 * Real.log l + C₀ :=
  ⟨0, fun x l hl _ => by
    have h := faltingsHeightJ_quotLCyclicJ_uncond x hl
    linarith⟩

/-- ★★★`α` の側（葉 3）を受ける。 -/
theorem imageContainsSL2J_torsionExt (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hl5 : 5 ≤ l) (hm : x.rep.toSSCurve.HasMultRed)
    (hp : x.rep.toSSCurve.PrimeToLocalHeights l)
    (hc : ¬ HasLCyclicJ x.rep.toSSCurve l) : ImageContainsSL2J x.rep.toSSCurve l := by
  sorry

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**葉 4 が閉じたので、これは `sorry` ではなくなった**（2026-08-31、第 784）。

`hram : l ∉ ramPrimes` は「`l` は `disc L` を割らない」を含む——
`ramPrimes` の定義に `(disc L).natAbs.primeFactors` が入っているからである。 -/
theorem imageSurjectiveJ_of_containsSL2' (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hram : x.PrimeToRamification l) (h : ImageContainsSL2J x.rep.toSSCurve l) :
    ImageSurjectiveJ x.rep.toSSCurve l := by
  haveI : Fact l.Prime := ⟨hl⟩
  have hunram : ¬ (l : ℤ) ∣ NumberField.discr x.rep.toSSCurve.fld := by
    intro hdvd
    refine hram (Finset.mem_union.2 (Or.inl (Finset.mem_union.2 (Or.inr ?_))))
    refine Nat.mem_primeFactors.2 ⟨hl, ?_, ?_⟩
    · have := Int.natAbs_dvd_natAbs.2 hdvd
      simpa using this
    · exact Int.natAbs_ne_zero.2 (NumberField.discr_ne_zero (K := x.rep.toSSCurve.fld))
  exact imageSurjectiveJ_of_cyclotomic x.rep.toSSCurve l
    (fun u => cyclotomic_det_surjective x.rep.toSSCurve l hunram u) h

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`EllModuliData` の witness**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★50 個近い欄のうち、`sorry` は上の 5 本だけである。 -/
noncomputable def ellModuliWitness : EllModuliData where
  EllClass := ℂ
  Curve := RealizedClass
  cls := RealizedClass.cls
  degOfDefinition := RealizedClass.degOfDefinition
  faltingsHeight := faltingsHeightJ
  HasPotMultRed := fun x => x.rep.toSSCurve.HasPotMultRed
  PrimeToLocalHeights := fun x l => x.rep.toSSCurve.PrimeToLocalHeights l
  CompactlyBounded := CompactlyBoundedJ
  GaloisFinite := GaloisFiniteJ
  ImageContainsSL2 := fun x l => ImageContainsSL2J x.rep.toSSCurve l
  degInf := degInfJ
  htInf := fun j => 12 * faltingsHeightJ j
  logDiffMell := logDiffMellJ
  degLe := degLeJ
  SemiStable := fun x => x.rep.toSSCurve.SemiStable
  HasLCyclic := fun x l => HasLCyclicJ x.rep.toSSCurve l
  MinimalField := fun x => MinimalFieldJ x.rep
  ImageSurjective := fun x l => ImageSurjectiveJ x.rep.toSSCurve l
  PrimeToRamification := RealizedClass.PrimeToRamification
  HasMultRed := fun x => x.rep.toSSCurve.HasMultRed
  PrimeToMultPrimes := fun x l => x.rep.PrimeToMultPrimes l
  degInf_le_htInf := degInfJ_sub_htInfJ_le
  htInf_bdeq_faltings := ⟨0, fun x => by simp⟩
  faltingsHeight_bddBelow := faltingsHeightJ_bddBelow
  northcott := fun C d _ => northcottJ C d
  quotLCyclic := quotLCyclicJ
  degInf_quotLCyclic := fun E l hl hcyc hpr => degInfJ_quotLCyclicJ E l hl hcyc hpr
  faltingsHeight_quotLCyclic := faltingsHeightJ_quotLCyclicJ
  degOfDefinition_pos := RealizedClass.degOfDefinition_pos
  primeToLocalHeights_of_lt := by
    intro E l hl _ hlt
    refine E.rep.toSSCurve.primeToLocalHeights_of_lt l hl ?_
    have h2 : degInfJ E.rep.toSSCurve.j = degInfJ E.cls := E.degInfJ_rep
    rw [h2]
    exact hlt
  lcyclicExc := lcyclicExcJ
  galoisFinite_lcyclicExc := galoisFiniteJ_lcyclicExcJ
  mem_lcyclicExc := by
    intro C eps KV E l hl hss hcyc hpr hor
    exact ⟨E, l, rfl, hl, hss, hcyc, hpr, hor⟩
  noMultRedExc := noMultRedExcJ
  galoisFinite_noMultRedExc := galoisFiniteJ_noMultRedExcJ
  mem_noMultRedExc := by
    intro KV E hKV hnm
    exact absurd E.rep.multRed hnm
  galoisFinite_union := galoisFiniteJ_union
  torsionExt := fun x => x
  cls_torsionExt := fun _ => rfl
  degOfDefinition_torsionExt := by
    intro E
    have := E.degOfDefinition_pos
    omega
  semiStable_torsionExt := fun E => E.rep.toSSCurve.semiStable_all
  hasMultRed_torsionExt := fun _ h => h
  primeToLocalHeights_torsionExt := fun _ _ h _ => h
  imageContainsSL2_of_torsionExt := fun E l hl hl5 hm hp hc =>
    imageContainsSL2J_torsionExt E l hl hl5 hm hp hc
  imageSurjective_of_containsSL2 := fun E l hl hram h =>
    imageSurjectiveJ_of_containsSL2' E l hl hram h
  compactlyBounded_empty := compactlyBoundedJ_empty
  multCard := fun x => x.rep.multCard
  multCard_pos := fun x => x.rep.multCard_pos
  multPrime := fun x j => x.rep.multPrime j
  multPrime_prime := fun x j => x.rep.multPrime_prime j
  localHt := fun x j => x.rep.localHt j
  localHt_pos := fun x j => x.rep.localHt_pos j
  sum_localHt_eq := RealizedClass.sum_localHt_eq
  badPrimes := fun x => x.rep.badPrimes
  badPrimes_prime := fun x => x.rep.badPrimes_prime
  sum_log_badPrimes_le := fun x => x.rep.sum_log_badPrimes_le
  primeTo_badPrimes := fun x l hl h => x.rep.primeTo_badPrimes l hl h
  ramPrimes := RealizedClass.ramPrimes
  ramPrimes_prime := RealizedClass.ramPrimes_prime
  badPrimes_subset_ramPrimes := RealizedClass.badPrimes_subset_ramPrimes
  primeTo_ramPrimes := RealizedClass.primeTo_ramPrimes
  sum_log_ramPrimes_le := by
    intro E
    have h := E.sum_log_ramPrimes_le
    rw [logDiffMellJ_eq]
    exact h

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★witness へ流し込む -/

/-- ★★★★★★★★★★★★★★★★★★★★**`Lemma 3.7` を witness で**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★これで `Lemma 3.7` が**具体的な楕円曲線の言葉**になった。
☆残る `sorry` は witness の 5 本だけである。 -/
theorem lemma_3_7_witness (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : RealizedClass) (l : ℕ), Nat.Prime l → E.rep.toSSCurve.SemiStable →
        ∀ (condA condB : Prop),
          (condA ↔ (100 * (E.degOfDefinition : ℝ)
                      * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ eps)
                        ≤ (l : ℝ) ∧ E.rep.toSSCurve.HasMultRed)) →
          (condB ↔ (E.cls ∈ KV ∧ E.rep.toSSCurve.PrimeToLocalHeights l)) →
          (condA → E.rep.toSSCurve.PrimeToLocalHeights l)
        ∧ (condB → E.cls ∉ Exc → E.rep.toSSCurve.HasMultRed)
        ∧ ((condA ∨ condB) → HasLCyclicJ E.rep.toSSCurve l → E.cls ∈ Exc) :=
  lemma_3_7 ellModuliWitness KV hKV eps heps

/-- ★★★★★★★★★★★★★★★★★★★★★★**`Theorem 3.8` を witness で**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem theorem_3_8_witness (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : RealizedClass) (l : ℕ), Nat.Prime l → E.cls ∉ Exc →
        ((23040 * 100 * (E.degOfDefinition : ℝ)
              * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ E.rep.toSSCurve.HasPotMultRed)
          ∨ (E.cls ∈ KV ∧ E.rep.toSSCurve.PrimeToLocalHeights l ∧ Nat.Coprime l 30)) →
        ImageContainsSL2J E.rep.toSSCurve l :=
  theorem_3_8 ellModuliWitness KV hKV ε hε

/-- ★★★★`Corollary 4.3` も witness で通ることの確認（型検査のみ）。 -/
example (eps : ℝ) (heps : 0 < eps) := cor_4_3 ellModuliWitness eps heps

/-- ★★★★`Corollary 4.4` も witness で通ることの確認（型検査のみ）。 -/
example (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) := cor_4_4 ellModuliWitness KV hKV

def lemma_3_7_witness.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(witness へ流し込んだ形——具体的な楕円曲線の言葉で)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_witness.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ellModuliWitness(残る sorry は 5 本の葉)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.ellModuliWitness") 1 ]

def theorem_3_8_witness.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(witness へ流し込んだ形——具体的な楕円曲線の言葉で)",
    sectionId := "genell-thm-3-8" }

def theorem_3_8_witness.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ellModuliWitness(残る sorry は 5 本の葉)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.ellModuliWitness") 1 ]

/-! ## ★出典の紐付け(`.src`) -/

def galoisFiniteJ_lcyclicExcJ_of_veluQuotOK.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(lcyclicExcJ は VeluQuotOK さえあれば Galois-有限)",
    sectionId := "genell-lemma-3-7" }

def galoisFiniteJ_lcyclicExcJ_of_veluQuotOK.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_galoisFinite_lcyclic(第 1239、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_galoisFinite_lcyclic") 1,
    .citation "[ABC3]" "galoisFiniteJ_subset(第 1238、無条件)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galoisFiniteJ_subset") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1341）**——界面の `galoisFinite_lcyclicExc` の " ++
       "`C` の量化を `∃ C₀, ∀ C ≥ C₀` に直した。" ++
       "☆かつては `∀ C` であり、`C` をいくらでも小さく取れるので" ++
       "**作れない主張**だった（条件 (a) が常に満たされてしまう）。" ++
       "★消費側（`Section3.lean` の `lemma_3_7`）は `C` を自分で選ぶので波及はない。") 3 ]

def lcyclicExcJ.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(lcyclicExc 欄——mem_lcyclicExc の条件そのものを集めた集合)",
    sectionId := "genell-lemma-3-7" }

def lcyclicExcJ.needs : List ProofObligation := []

def galoisFiniteJ_lcyclicExcJ.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(lcyclicExc は Galois-finite——Lemma 3.5 の (†) を受ける)",
    sectionId := "genell-lemma-3-7" }

def galoisFiniteJ_lcyclicExcJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Lemma 3.5 の (†)((l/14)·deg∞ ≤ ht^Falt + 2log l + C′)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hdag_of_velu") 3,
    .implicitStep
      ("★(a)(b) 両側の高さ評価は済んでいる(htFalt_le_of_condA・htFalt_le_of_condB)。" ++
       "有限性は northcottJ(§9-1179、第 753、★無条件)から出る") 4 ]

def degInfJ_quotLCyclicJ_ge_of_exists.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(商の類が存在すれば l·deg∞(E) ≤ deg∞(E′)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ_ge_of_exists.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1,
    .citation "[ABC3]" "degInfOf_ge_of_local / minDeltaExp_le_of_bad_delta(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_ge_of_local") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1340）**——界面の `degInf_quotLCyclic` を" ++
       "**等式から不等式に弱めた**ことで、" ++
       "良い素点の半安定性（同種で良還元が保たれること）が不要になった。" ++
       "☆消費側（`Section3.lean` の `lemma_3_5`）が使うのはこの向きだけである。") 3 ]

def degInfJ_quotLCyclicJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(degInf_quotLCyclic 欄——deg∞(E′) = l·deg∞(E))",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_large(第 1044、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_large") 1,
    .citation "[ABC3]" "degInfOf_eq_of_local(局所の関係を足し上げる、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_eq_of_local") 1,
    .implicitStep
      ("★★**2026-09-01（第 1086）の測定**——" ++
       "本欄の数学はすでにある。" ++
       "☆`IsMuAtBadPrimes`（第 1044、`sorry`-free）が悪い素点で " ++
       "`Δ_min(E′) = l·Δ_min(E)` を与え、半安定なら良い素点では両辺 `0`、" ++
       "`degInfOf_eq_of_local` がそれを足し上げる。" ++
       "★残っているのは**実現の橋**だけである——" ++
       "`quotLCyclicJ x l` と `veluQuotientFull` を同定し、" ++
       "`degInfJ` と `degInfOf` を同定する段。") 6 ]

def faltingsHeightJ_quotLCyclicJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(faltingsHeight_quotLCyclic 欄——ht^Falt(E′) ≤ ht^Falt(E) + 2log l + C₀)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quotLCyclicJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_veluQuotientFull_le(§9-1160、第 704、★無条件)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le") 3,
    .implicitStep "☆quotLCyclicJ の定義と Vélu の商を突き合わせる段" 5 ]

def imageContainsSL2J_torsionExt.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(imageContainsSL2_of_torsionExt 欄)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_torsionExt.needs : List ProofObligation :=
  [ .citation "[ABC3]" "alpha_in_modl_image(Skeleton/GenEll/GaloisLocal.lean、局所理論の行列表示)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.alpha_in_modl_image") 10,
    .implicitStep
      ("★群論(Lemma 3.1, (iv))と位相(連続性・閉性)は済んでいる" ++
       "——imageContainsSL2J_of_alpha'(§9-1200、第 774)") 1 ]

def imageSurjectiveJ_of_containsSL2'.src : Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(imageSurjective_of_containsSL2 欄)",
    sectionId := "genell-cor-4-3" }

def imageSurjectiveJ_of_containsSL2'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cyclotomic_det_surjective(Skeleton/GenEll/GaloisLocal.lean、円分指標の全射性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.cyclotomic_det_surjective") 1,
    .implicitStep
      ("★群論と逆極限の段は済んでいる" ++
       "——imageSurjectiveJ_of_cyclotomic(第 765)と " ++
       "cyclotomicCharacter_surjective_of_mod(§9-1199、第 773)") 1 ]

def ellModuliWitness.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の witness——sorry は 5 本の葉だけ)",
    sectionId := "genell-prop-3-4" }

end ABC3.Skeleton.GenEll
