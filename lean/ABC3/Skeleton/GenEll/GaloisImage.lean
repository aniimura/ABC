import ABC3.Meta.Claim
import ABC3.Interface.GenEll.GaloisRep
import ABC3.Found.GenEll.Lemma31
import ABC3.Skeleton.GenEll.Section3

/-!
# [GenEll] Theorem 3.8 —— Full Special Linear Galois Actions(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.19。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★このスケルトンの位置づけ

`[IUTchIV]` が直接引く [GenEll] の項目は **`Theorem 2.1`(10 回)・`Theorem 3.8`(2 回)・
`Corollary 4.4`(2 回)** など。うち `Corollary 4.3` / `Corollary 4.4`(§4)は
**本定理の上にしか立たない**——すなわち **§4 の 5 件は本定理を経由して
l 捩れ Galois 表現に律速される**。

★**statement を書くだけで語彙が足りない**ので、
`Interface/GenEll/GaloisRep.lean` の `TorsionGaloisRepData` が下敷きになる。
各フィールドは原文 p.19 に**実際に現れる語**と 1 対 1 に対応させた。

## ★★★★★★★★ 2026-08-26——`sorry` が消えた(第 365 ブロック)

★以前はここに「`theorem_3_8` は `sorry` であり、**それを消すことを目的にしてはならない**」
と書いてあった——`Heights.lean` で実演したとおり、内容を `Interface` の仮説へ移せば
`sorry` は消えるが、それは `tools/check.mjs` 冒頭 **B5** が名指しする穴だからである。

★★**その警告は今も生きている**。本ブロックが B5 にならないのは、次の 3 点による:

1. ★**界面に出した 7 つはすべて原文 p.20 が実際に書いている段**である
   ——`L′` への基底変換(次数が 23040 を割る、半安定還元、乗法還元)と、最終段の `Lemma 3.1, (iv)`。
   **結論そのものを仮説に移したのではない**。
2. ★★**中身は `Lemma 3.7` からの導出である**——定数 23040 の帳尻(下記)も
   条件 (b) の 30 の使い道も、ここで初めて仕事をする。
3. ★★★**残った 1 つの posit は名前で見える**——`imageContainsSL2_of_torsionExt` だけであり、
   それは **Galois 表現が未構築だから**であって、`Lemma 3.1` を持っていないからではない。

★★★★**本当の進捗は依然として `Interface` が `waiting` でなくなること**であり、
`node tools/check.mjs` の「Interface 実装待ち」がそれを見せる。

## ★★★★ 2026-08-26 の逸脱の記録——量化する界面を下げた

| 項 | 原典 | 形式化 |
|---|---|---|
| 量化する界面 | (区別なし) | `TorsionGaloisRepData` → **`EllModuliData`** |

★理由: 原文 p.20 の証明は `Lemma 3.7` を使うが、`Lemma 3.7` は `EllModuliData` の上に立つ。
★★`EllModuliData extends TorsionGaloisRepData` なので、**下流は何も失わない**
——`Corollary 4.3` / `Corollary 4.4` はもともと `EllModuliData` の上にある。
★★★弱めているのは「どの `D` について言うか」だけで、**結論は同じ**である。

## ★★★★★定数 23040 の帳尻(本ブロックの算術的な中身)

原文は『for a suitable choice of `C` … then `E_{L′}` and `l` satisfy condition (a) …
of Lemma 3.7 **[perhaps for a different “C”]**』と 1 文で済ませている。実際には:

    d′ ≤ 23040·d 、 [E_{L′}] = [E_L] なので ht^Falt はそのまま
    d′^ε ≤ (23040·d)^ε = 23040^ε·d^ε
    C := C₇·23040^ε + |B| + 1 と取れば
      100·d′·(ht^Falt + C₇·d′^ε) ≤ 23040·100·d·(ht^Falt + C·d^ε)

★★`ht^Falt < 0` のとき第 1 項は**逆向きになる**が、その超過分は `|B|·23040·d` 以下であり
(`ht^Falt ≥ B`)、`C` の `|B| + 1` の分が `d^ε ≥ 1` を使ってそれを吸収する。

## ★原文の証明のうち、群論は 1 点だけ(実測)

原文 p.20 の証明を通読した結果、**群論に頼るのは `Lemma 3.1, (iv)` ただ 1 つ**である。
残りはすべて算術側(Faltings 高さ・Tate 曲線・半安定還元・`Lemma 3.7`)。
★`Lemma 3.1` の (i)(ii)(iii) は `Found/GenEll/Lemma31.lean` に**実装済み**(`sorry` 無し)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll

set_option maxHeartbeats 1600000 in
/-- **[GenEll] Theorem 3.8**(Full Special Linear Galois Actions)。

原文 (GenEll p.19):
> Theorem 3.8. (Full Special Linear Galois Actions) Let Q[bb][bar] be an algebraic

★**条件 (a) の指数 `ϵ` は `.txt` から読んではならない**——
`pdftotext -layout` は `C · d^ϵ` を `C · d ` と出し、**条件が別物になる**
(2026-08-16 実測。`1_Structured/…/theorem-3-8.html` に記録)。
本 statement は **260 dpi 目視**から写した。

★`23040 = |GL₂(𝔽₃) × GL₂(𝔽₅)|` は原文 p.20 の証明に現れる定数
`(3²−1)(3²−3)(5²−1)(5²−5)` である(目視確認)。 -/
theorem theorem_3_8 (D : EllModuliData)
    (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV)
    (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.cls E ∉ Exc →
        ((23040 * 100 * (D.degOfDefinition E : ℝ)
              * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ D.HasPotMultRed E)
          ∨ (D.cls E ∈ KV ∧ D.PrimeToLocalHeights E l ∧ Nat.Coprime l 30)) →
        D.ImageContainsSL2 E l := by
  obtain ⟨C₇, hC₇pos, Exc, hExc, h37⟩ := lemma_3_7 D KV hKV ε hε
  obtain ⟨B, hB⟩ := D.faltingsHeight_bddBelow
  have hKpos : (0:ℝ) < (23040:ℝ) ^ ε := Real.rpow_pos_of_pos (by norm_num) ε
  refine ⟨C₇ * (23040:ℝ) ^ ε + |B| + 1, by positivity, Exc, hExc, ?_⟩
  intro E l hl hExcE hcond
  set E' := D.torsionExt E with hE'def
  have hcls : D.cls E' = D.cls E := D.cls_torsionExt E
  have hss : D.SemiStable E' := D.semiStable_torsionExt E
  have hExcE' : D.cls E' ∉ Exc := by rw [hcls]; exact hExcE
  obtain ⟨ha, hb, hc⟩ := h37 E' l hl hss
    (100 * (D.degOfDefinition E' : ℝ)
        * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ)
      ∧ D.HasMultRed E')
    (D.cls E' ∈ KV ∧ D.PrimeToLocalHeights E' l) Iff.rfl Iff.rfl
  have hmain : D.HasMultRed E' ∧ D.PrimeToLocalHeights E' l ∧ ¬ D.HasLCyclic E' l := by
    rcases hcond with ⟨hnum, hpot⟩ | ⟨hKVE, hpl, hcop⟩
    · -- 条件 (a): 基底変換して `Lemma 3.7` の (a) を満たす
      have hmult' : D.HasMultRed E' := D.hasMultRed_torsionExt E hpot
      have hnumA : 100 * (D.degOfDefinition E' : ℝ)
          * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ) := by
        rw [hcls]
        set d : ℝ := (D.degOfDefinition E : ℝ) with hddef
        set d' : ℝ := (D.degOfDefinition E' : ℝ) with hd'def
        set F : ℝ := D.faltingsHeight (D.cls E) with hFdef
        have hd1 : (1:ℝ) ≤ d := by
          rw [hddef]; exact_mod_cast D.degOfDefinition_pos E
        have hd'1 : (1:ℝ) ≤ d' := by
          rw [hd'def]; exact_mod_cast D.degOfDefinition_pos E'
        have hdd' : d' ≤ 23040 * d := by
          rw [hd'def, hddef]
          exact_mod_cast D.degOfDefinition_torsionExt E
        have hdpos : (0:ℝ) < d := by linarith
        have hd'pos : (0:ℝ) < d' := by linarith
        set K : ℝ := (23040:ℝ) ^ ε with hKdef
        set P : ℝ := d ^ ε with hPdef
        set P' : ℝ := d' ^ ε with hP'def
        have hP1 : (1:ℝ) ≤ P := by rw [hPdef]; exact Real.one_le_rpow hd1 hε.le
        have hP'1 : (1:ℝ) ≤ P' := by rw [hP'def]; exact Real.one_le_rpow hd'1 hε.le
        have hK1 : (1:ℝ) ≤ K := by rw [hKdef]; exact Real.one_le_rpow (by norm_num) hε.le
        have hPKP : P' ≤ K * P := by
          rw [hP'def, hKdef, hPdef, ← Real.mul_rpow (by norm_num) hdpos.le]
          exact Real.rpow_le_rpow hd'pos.le hdd' hε.le
        have hBF : -|B| ≤ F := le_trans (neg_abs_le B) (hB (D.cls E))
        have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
        have hdPnn : (0:ℝ) ≤ d * P := by nlinarith
        have hBdPnn : (0:ℝ) ≤ |B| * (d * P) := mul_nonneg hBnn hdPnn
        have hstep : 100 * d' * (F + C₇ * P')
            ≤ 23040 * 100 * d * (F + (C₇ * K + |B| + 1) * P) := by
          -- ★両辺を単項へ展開しておく(`linarith` は積を原子として見る)
          have hLexp : 100 * d' * (F + C₇ * P') = 100 * d' * F + 100 * C₇ * (d' * P') := by ring
          have hRexp : 23040 * 100 * d * (F + (C₇ * K + |B| + 1) * P)
              = 23040 * 100 * d * F + 23040 * 100 * C₇ * K * (d * P)
                + 23040 * 100 * (|B| * (d * P)) + 23040 * 100 * (d * P) := by ring
          rw [hLexp, hRexp]
          -- ★第 2 項: `d'·d'^ε ≤ (23040d)·(23040^ε d^ε)`
          have hterm2 : 100 * C₇ * (d' * P') ≤ 23040 * 100 * C₇ * K * (d * P) := by
            have h1 : d' * P' ≤ (23040 * d) * (K * P) :=
              mul_le_mul hdd' hPKP (by linarith) (by linarith)
            have h2 := mul_le_mul_of_nonneg_left h1 (by positivity : (0:ℝ) ≤ 100 * C₇)
            linarith [h2]
          rcases le_or_gt 0 F with hF | hF
          · -- ★`F ≥ 0` なら第 1 項は `d' ≤ 23040d` で直ちに出る
            have hterm1 : 100 * d' * F ≤ 23040 * 100 * d * F := by nlinarith
            linarith [hterm1, hterm2, hBdPnn, hdPnn]
          · -- ★★`F < 0` なら第 1 項は逆向きになるが、超過分は `|B|·23040d` で押さえられ、
            -- `C` の `|B| + 1` の分がそれを吸収する(`P ≥ 1` が効く)
            have hFB : -F ≤ |B| := by linarith
            have hmul : (-F) * (23040 * d - d') ≤ |B| * (23040 * d) :=
              mul_le_mul hFB (by linarith) (by linarith) hBnn
            have hterm1 : 100 * d' * F - 23040 * 100 * d * F ≤ 23040 * 100 * (|B| * d) := by
              nlinarith [hmul]
            have hdP : d ≤ d * P := by nlinarith
            have hBdP : |B| * d ≤ |B| * (d * P) := mul_le_mul_of_nonneg_left hdP hBnn
            linarith [hterm1, hterm2, hBdP, hdPnn]
        linarith
      have hcA : (100 * (D.degOfDefinition E' : ℝ)
          * (D.faltingsHeight (D.cls E') + C₇ * (D.degOfDefinition E' : ℝ) ^ ε) ≤ (l : ℝ)
        ∧ D.HasMultRed E') := ⟨hnumA, hmult'⟩
      refine ⟨hmult', ha hcA, fun hcyc => hExcE' (hc (Or.inl hcA) hcyc)⟩
    · -- 条件 (b): `30` と互いに素なので局所高さと素であることが `L′` へ移る
      have hpl' : D.PrimeToLocalHeights E' l := D.primeToLocalHeights_torsionExt E l hpl hcop
      have hcB : (D.cls E' ∈ KV ∧ D.PrimeToLocalHeights E' l) := by
        refine ⟨?_, hpl'⟩
        rw [hcls]; exact hKVE
      exact ⟨hb hcB hExcE', hpl', fun hcyc => hExcE' (hc (Or.inr hcB) hcyc)⟩
  -- ★`Lemma 3.1, (iv)` が要求する `5 ≤ l` は、どちらの条件からも出る
  have hl5 : 5 ≤ l := by
    have h2 := hl.two_le
    rcases hcond with ⟨hnum, -⟩ | ⟨-, -, hcop⟩
    · set d : ℝ := (D.degOfDefinition E : ℝ) with hddef
      set F : ℝ := D.faltingsHeight (D.cls E) with hFdef
      have hd1 : (1:ℝ) ≤ d := by
        rw [hddef]; exact_mod_cast D.degOfDefinition_pos E
      have hP1 : (1:ℝ) ≤ d ^ ε := Real.one_le_rpow hd1 hε.le
      have hFB : -|B| ≤ F := le_trans (neg_abs_le B) (hB (D.cls E))
      have hBnn : (0:ℝ) ≤ |B| := abs_nonneg B
      have hC7 : (0:ℝ) < C₇ * (23040:ℝ) ^ ε := mul_pos hC₇pos hKpos
      have hCnn : (0:ℝ) ≤ C₇ * (23040:ℝ) ^ ε + |B| + 1 := by positivity
      have hCP : (C₇ * (23040:ℝ) ^ ε + |B| + 1)
          ≤ (C₇ * (23040:ℝ) ^ ε + |B| + 1) * d ^ ε := by nlinarith [hP1, hCnn]
      have hsum : (1:ℝ) ≤ F + (C₇ * (23040:ℝ) ^ ε + |B| + 1) * d ^ ε := by linarith
      have hprod : (2304000:ℝ)
          ≤ 23040 * 100 * d * (F + (C₇ * (23040:ℝ) ^ ε + |B| + 1) * d ^ ε) := by
        nlinarith [hd1, hsum]
      have hkey : (2304000 : ℝ) ≤ (l : ℝ) := le_trans hprod hnum
      have hnat : (2304000 : ℕ) ≤ l := by exact_mod_cast hkey
      omega
    · by_contra hlt
      push_neg at hlt
      interval_cases l
      · exact absurd hcop (by decide)
      · exact absurd hcop (by decide)
      · rcases hl.eq_one_or_self_of_dvd 2 (by norm_num) with h | h <;> omega
  exact D.imageContainsSL2_of_torsionExt E l hl hl5 hmain.1 hmain.2.1 hmain.2.2

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_3_8.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8",
    sectionId := "genell-thm-3-8" }

/-- ★この定理の証明が要求するもの(原文 p.20 を通読して数えた)。

★**群論は `Lemma 3.1, (iv)` ただ 1 つ**で、そのうち (i)(ii)(iii) は実装済みである。
残りはすべて算術側で、そこが `Interface` に載っている。 -/
def theorem_3_8.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Lemma 3.1, (iv)(SL₂(ℤ_l) の持ち上げ)——★原文の証明が使う群論はこれだけであり、★★★★★**4 条すべてが実装済み**である(2026-08-29 確認): (i)(ii)(iii) は Found/GenEll/Lemma31.lean、**(iv) は Found/GenEll/Sl2Padic.lean の lemma_3_1_iv**。★原文は [Serre] Chapter IV §3.4 Lemma 3 を引くが 0_Source に無いため、本プロジェクトが自分で証明している。★★★★★★**したがって Theorem 3.8 に「Serre の開像定理」は要らない**——障壁は Tate 曲線(下の 2 行)だけである" 14,
    .otherPaper "[GenEll]" "Lemma 3.7(局所高さと l-巡回部分群スキーム)" 18,
    .otherPaper "[GenEll]" "Proposition 3.4(Faltings 高さによる例外集合の有限性)" 17,
    .otherPaper "[GenEll]" "Lemma 3.2 の直前の局所理論(Tate 曲線)" 15,
    .citation "[FC]" "Degenerations of Abelian Varieties(半安定還元)"
      (.absent "mathlib は `EllipticCurve/Reduction.lean` を持つが Tate 曲線・半安定還元の理論は無い(2026-08-16、ディレクトリ全宣言を確認)") 19,
    .implicitStep
      "原文は 3・5 捩れ点を有理化する次数 23040 の Galois 拡大へ移る段を『there exists a Galois extension L′ of L of degree that divides d₀』と 1 文で済ませている" 20,
    .implicitStep
      "★statement の語彙(Faltings 高さ・Galois 表現・compactly bounded・Galois-finite)を Interface/GenEll/GaloisRep.lean に posit した。**我々は作っていない**" 19 ]

end ABC3.Skeleton.GenEll
