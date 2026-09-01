/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.Lemma37A
import ABC3.Found.GaloisRep.Lemma37CFull
import ABC3.Found.GaloisRep.Lemma37CondBFull
import ABC3.Meta.Claim

/-!
# 第 1151 ブロック —— **`Lemma 3.7` を `Found` で束ねる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★これは何か

`Skeleton/GenEll/Section3.lean` の `lemma_3_7` は抽象界面 `EllModuliData D` の上の主張で、
`D` の欄（`primeToLocalHeights_of_lt`・`lcyclicExc`・`northcott` 等）を**外から受けている**。

★本ファイルは同じ主張を **`Found` の具体的な語彙だけ**で述べ、証明する。
☆`EllModuliData` の witness をすべて組み上げる必要は無い——`Lemma 3.7` が使う欄は
すでに `Found/` に揃っているからである。

| 主張 | 中身 | 出どころ |
|---|---|---|
| 第 1 | 条件 (a) なら `l` は局所高さと素 | `lemma_3_7_a_coprime`（★無条件） |
| 第 2 | 条件 (b) かつ `∉ Exc` なら乗法還元 | ★`DegCurve` は定義から乗法還元を持つので自明 |
| 第 3 | (a)∨(b) かつ `l`-巡回なら `∈ Exc` | `htFalt_le_of_condA_lcyclic`（第 1088）＋ `htFalt_le_of_condB_lcyclic`（第 1150） |
| `Exc` の有限性 | `ht^Falt` が有界な類は各次数で有限 | `galoisFiniteJ_htFalt_le`（第 756）＝ `northcottJ`（★無条件） |

## ★★★★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 影響 |
|---|---|---|---|
| 曲線の族 | 任意の楕円曲線 | ☆`DegCurve`（**半安定**かつ乗法還元を持つもの） | ★原文は `L(E[3],E[5])` への基底変換で一般の場合を半安定へ帰着させる。その段を後回しにしている |
| `H_L ⊆ E_L` は `l`-巡回部分群スキーム | スキーム | ☆`HasLCyclicVelu`——**`L` 有理点 `Q` で `addOrderOf Q = l`**、その Vélu の商が楕円曲線で半安定 | ★階数 1 の部分群を生成元で書いたもの。☆商の半安定性は原文が「同種なので自動」と括弧で述べる |
| `Galois-finite` | 有限個の Galois 軌道の合併 | ☆`GaloisFiniteJ`——各次数 `d` で `S ∩ M_ell(ℚ̄)^{≤d}` が有限 | ★原文より**弱い**が、`Exc` は結論側に現れるので安全である |
| `compactly bounded` | 有限個の素点 `V` とコンパクト `K_v` | ☆`CompactlyBoundedJ`——`h_∞(j)` が有界 | ★アルキメデス素点での内容そのもの |

☆これらはすべて `Found/GenEll/EllModuliObjects.lean` で先に決めた設計であり、
本ファイルで新たに増やした逸脱は無い。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Meta

/-! ## ★★★★★★局所高さと素なら `l ∤ jExp` -/

/-- ★★★★★★★★**半安定なら「局所高さと素」は「`l ∤ v_p(j)`」と同じ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆半安定なら `v_p(Δ_min) = max(0, −v_p(j))` なので、`v_p(j) < 0` の素点では
局所高さは `−v_p(j)` そのものである。★`Lemma 3.5` が受ける形に翻訳する橋である。 -/
theorem SSCurve.not_dvd_jExp_of_primeToLocalHeights (E : SSCurve) {l : ℕ} (hl : Nat.Prime l)
    (h : E.PrimeToLocalHeights l) :
    ∀ p : HeightOneSpectrum (𝓞 E.fld), jExp p E.W < 0 → ¬ ((l : ℤ) ∣ jExp p E.W) := by
  intro p hneg hdvd
  have hmd : minDeltaExp p E.W = max 0 (-jExp p E.W) :=
    minDeltaExp_eq_maxJ_of_semistable p E.W (E.ss p)
  have hmax : max (0 : ℤ) (-jExp p E.W) = -jExp p E.W := max_eq_right (by omega)
  rw [hmax] at hmd
  have hne : E.localHeightAt p ≠ 0 := by
    show minDeltaExp p E.W ≠ 0
    rw [hmd]; omega
  have hcop := h p hne
  have hnn : 0 ≤ minDeltaExp p E.W := minDeltaExp_nonneg p E.W
  have hdvdZ : (l : ℤ) ∣ minDeltaExp p E.W := by rw [hmd]; exact dvd_neg.2 hdvd
  have hcast : ((minDeltaExp p E.W).toNat : ℤ) = minDeltaExp p E.W := Int.toNat_of_nonneg hnn
  have hdvdN : l ∣ (minDeltaExp p E.W).toNat := by
    have : ((l : ℕ) : ℤ) ∣ (((minDeltaExp p E.W).toNat : ℕ) : ℤ) := by rw [hcast]; exact hdvdZ
    exact_mod_cast this
  rw [Nat.Prime.coprime_iff_not_dvd hl] at hcop
  exact hcop hdvdN

/-! ## ★★★★★★★★`l`-巡回部分群を Vélu の商で書いたもの -/

/-- ★★★★★★★★★★**`HasLCyclic` 欄を `Lemma 3.5` が受ける形で**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

☆原文『`E_L` admits an `l`-cyclic subgroup scheme `H_L ⊆ E_L`』を
**生成元 `Q`（`L` 有理点、位数 `l`）と、その Vélu の商**で書いたものである。
★商が楕円曲線であること・半安定であることを一緒に持たせるのは、
原文が「同種なので自動」と括弧で述べる段を仮説の形にしているためである
（`Found/GaloisRep/Lemma35Unconditional.lean` の逸脱表と同じ）。 -/
def HasLCyclicVelu (E : SSCurve) (l : ℕ) : Prop :=
  -- ★判定のインスタンスを `Classical` に固定する——`Lemma 3.5` 側の命題が
  -- 一般の体 `L` で `Classical.propDecidable` を拾っているので、
  -- ここで `Subtype` の判定を拾うと `Point` の群構造が別物になる。
  letI : DecidableEq E.fld := fun a b => Classical.propDecidable (a = b)
  ∃ Q : E.W.toAffine.Point, addOrderOf Q = l ∧
    (veluQuotientFull E.W (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).IsElliptic ∧
    ∀ p : HeightOneSpectrum (𝓞 E.fld),
      SemistableAt p (veluQuotientFull E.W (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q))))

/-! ## ★★★★★★★★★★★★★★★★★★★★`Lemma 3.7` そのもの -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Lemma 3.7**（Finite Exceptional Sets）——`Found` の語彙で。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

`K_V` を compactly bounded、`ϵ > 0` とすると、**定数 `C > 0` と
各次数で有限な例外集合 `Exc` が存在して**次を満たす:
`E` が半安定で乗法還元を持ち、`d ≝ [L:ℚ]`、`l` 素数として、

* (a) `l ≥ 100d·(ht^Falt + C·d^ϵ)` かつ `E` が乗法還元の素点を持つ
* (b) `[E] ∈ K_V` かつ `l` が乗法還元の全素点での局所高さと素

のいずれかを仮定すると:

1. (a) なら `l` は局所高さと**素**である。
2. (b) かつ `[E] ∉ Exc` なら `E` は乗法還元の素点を 1 つ以上持つ。
3. (a) か (b) が成り立ち、さらに `E` が `l`-巡回部分群を持つなら **`[E] ∈ Exc`**。

★★★★**2026-09-01（第 1151）**——外部からの入力は一切無い。
☆第 1 の主張は `lemma_3_7_a_coprime`（★無条件）、
第 2 の主張は `DegCurve` が定義から乗法還元を持つので自明、
第 3 の主張は `htFalt_le_of_condA_lcyclic`（第 1088）と
`htFalt_le_of_condB_lcyclic`（第 1150）である。
★`Exc` の有限性は `galoisFiniteJ_htFalt_le`（＝ `northcottJ`、★無条件）。 -/
theorem lemma_3_7 (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : DegCurve) (l : ℕ), Nat.Prime l →
        ∀ condA condB : Prop,
          (condA ↔ (100 * (E.deg : ℝ)
                      * (faltingsHeightJ E.j + C * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)
                    ∧ E.toSSCurve.HasMultRed)) →
          (condB ↔ (E.j ∈ KV ∧ E.toSSCurve.PrimeToLocalHeights l)) →
          (condA → E.toSSCurve.PrimeToLocalHeights l)
        ∧ (condB → E.j ∉ Exc → E.toSSCurve.HasMultRed)
        ∧ ((condA ∨ condB) → HasLCyclicVelu E.toSSCurve l → E.j ∈ Exc) := by
  obtain ⟨M, hM⟩ := hKV
  obtain ⟨Ca, hCa0, hCa⟩ := lemma_3_7_a_coprime eps heps
  obtain ⟨Cc, C'a, hCc0, hCc⟩ := htFalt_le_of_condA_lcyclic eps heps
  obtain ⟨C₂, C'b, _, hCb⟩ := htFalt_le_of_condB_lcyclic
  refine ⟨max Ca Cc, lt_of_lt_of_le hCa0 (le_max_left _ _),
    {x : ℂ | faltingsHeightJ x ≤ max C'a (|M + C₂| / 5 + 28 / 5 + 1.4 * C'b)},
    galoisFiniteJ_htFalt_le _, ?_⟩
  intro E l hl condA condB hcA hcB
  -- ☆語彙の橋
  have hFeq : faltingsHeightJ E.j = htFaltOf E.toSSCurve.fld E.toSSCurve.W :=
    faltingsHeightJ_eq E.toSSCurve
  have hdeq : (E.deg : ℝ) = (Module.finrank ℚ E.toSSCurve.fld : ℝ) := rfl
  have hd1 : (1 : ℝ) ≤ (E.deg : ℝ) := by exact_mod_cast E.deg_pos
  have hp0 : (0 : ℝ) ≤ (E.deg : ℝ) ^ eps := Real.rpow_nonneg (by linarith) eps
  -- ☆`C` を小さいものに置き換えてよい
  have hweak : ∀ C₀ : ℝ, C₀ ≤ max Ca Cc →
      (100 * (E.deg : ℝ) * (faltingsHeightJ E.j + max Ca Cc * (E.deg : ℝ) ^ eps) ≤ (l : ℝ)) →
      100 * (Module.finrank ℚ E.toSSCurve.fld : ℝ)
        * (htFaltOf E.toSSCurve.fld E.toSSCurve.W
            + C₀ * (Module.finrank ℚ E.toSSCurve.fld : ℝ) ^ eps) ≤ (l : ℝ) := by
    intro C₀ hle hA
    rw [← hdeq, ← hFeq]
    refine le_trans ?_ hA
    have hmul : C₀ * (E.deg : ℝ) ^ eps ≤ max Ca Cc * (E.deg : ℝ) ^ eps :=
      mul_le_mul_of_nonneg_right hle hp0
    have h100 : (0 : ℝ) ≤ 100 * (E.deg : ℝ) := by linarith
    exact mul_le_mul_of_nonneg_left (by linarith) h100
  -- ★第 1 の主張
  have hA1 : condA → E.toSSCurve.PrimeToLocalHeights l := by
    intro hc
    rw [hcA] at hc
    intro p hp
    exact hCa E.toSSCurve.fld E.toSSCurve.W l p hl (hweak Ca (le_max_left _ _) hc.1) hp
  refine ⟨hA1, fun _ _ => E.multRed, ?_⟩
  -- ★第 3 の主張
  letI : DecidableEq E.toSSCurve.fld := fun a b => Classical.propDecidable (a = b)
  intro hAB hcyc
  obtain ⟨Q, hQ, hell, hss'⟩ := hcyc
  haveI : (veluQuotientFull E.toSSCurve.W (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsElliptic := hell
  show faltingsHeightJ E.j ≤ max C'a (|M + C₂| / 5 + 28 / 5 + 1.4 * C'b)
  rw [hFeq]
  rcases hAB with hc | hc
  · -- ☆条件 (a)
    have hcA' := hcA.1 hc
    obtain ⟨p, hp⟩ := hcA'.2
    have hcop := E.toSSCurve.not_dvd_jExp_of_primeToLocalHeights hl (hA1 hc)
    have hmain := hCc E.toSSCurve.fld E.toSSCurve.W
      (veluQuotientFull E.toSSCurve.W (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))) l p hl hp
      (hweak Cc (le_max_right _ _) hcA'.1) Q hQ rfl E.toSSCurve.ss hss' hcop
    exact le_trans hmain (le_max_left _ _)
  · -- ☆条件 (b)
    have hcB' := hcB.1 hc
    have hcop := E.toSSCurve.not_dvd_jExp_of_primeToLocalHeights hl hcB'.2
    have harch : htArchJ E.toSSCurve.fld E.toSSCurve.W ≤ M := hM E.toSSCurve hcB'.1
    have hmain := hCb M E.toSSCurve.fld E.toSSCurve.W
      (veluQuotientFull E.toSSCurve.W (((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))) l hl E.toSSCurve.ss hss' harch Q hQ rfl hcop
    exact le_trans hmain (le_max_right _ _)

/-! ## ★出典の紐付け(`.src`) -/

def SSCurve.not_dvd_jExp_of_primeToLocalHeights.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(半安定なら「局所高さと素」は「l ∤ v_p(j)」と同じ。★無条件)",
    sectionId := "genell-lemma-3-7" }

def HasLCyclicVelu.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(l-巡回部分群を生成元と Vélu の商で書いたもの)",
    sectionId := "genell-lemma-3-7" }

def HasLCyclicVelu.needs : List ProofObligation :=
  [ .implicitStep
      ("★★**逸脱の記録**——原文の `H_L ⊆ E_L` は `l`-巡回**部分群スキーム**である。" ++
       "☆本定義は生成元を `L` 有理点に取り、さらに Vélu の商が楕円曲線で半安定であることを" ++
       "一緒に持たせている。★後者は原文が「同種なので自動」と括弧で述べる段であり、" ++
       "「同種なら半安定」の形式化は未了である。") 5 ]

def lemma_3_7.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7.needs : List ProofObligation :=
  [ .citation "[ABC3]" "lemma_3_7_a_coprime(第 1 の主張。★無条件、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.lemma_3_7_a_coprime") 1,
    .citation "[ABC3]" "htFalt_le_of_condA_lcyclic(第 3 の主張・条件 (a)、第 1088、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_le_of_condA_lcyclic") 1,
    .citation "[ABC3]" "htFalt_le_of_condB_lcyclic(第 3 の主張・条件 (b)、第 1150、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_le_of_condB_lcyclic") 1,
    .citation "[ABC3]" "galoisFiniteJ_htFalt_le(Exc の有限性＝ northcottJ。★無条件、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galoisFiniteJ_htFalt_le") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1151）の逸脱の記録**——曲線の族を `DegCurve`" ++
       "（半安定かつ乗法還元を持つもの）に制限している。" ++
       "☆原文は `L(E[3],E[5])` への基底変換で一般の場合を半安定へ帰着させるが、" ++
       "その段（`torsionExt` 欄）は後回しにしている。" ++
       "★`GaloisFinite` は「各次数で有限」と原文より弱く取っているが、" ++
       "`Exc` は結論側に現れるので安全である。") 8 ]

end ABC3.Found.GenEll
