#!/usr/bin/env python3
"""
Добавляет столбец pValue к merged_ligands_docking_*.csv файлам.
pValue = 9 - log10(standard_value), где standard_value в нМ.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np


def assign_activity_class(standard_value: float) -> str | None:
    """
    Классифицирует соединение по стандартному значению активности (Ki в нМ).
    
    Пороги (по просьбе пользователя):
      - standard_value <= 100 нМ       -> 'high'        (высокоактивное)
      - 100 < standard_value < 1000 нМ -> 'medium'      (среднеактивное)
      - standard_value >= 1000 нМ      -> 'low'         (неактивное)
    """
    if standard_value is None or np.isnan(standard_value):
        return None
    try:
        v = float(standard_value)
    except Exception:
        return None
    if v <= 100.0:
        return "high"
    if v < 1000.0:
        return "medium"
    return "low"

def add_pvalue_column(input_dir: Path, output_dir: Path):
    """
    Добавляет столбец pValue к CSV файлам.
    
    Args:
        input_dir: директория с исходными merged_ligands_docking_*.csv файлами
        output_dir: директория для сохранения файлов с pValue
    """
    print("=" * 80)
    print("ДОБАВЛЕНИЕ СТОЛБЦА pValue")
    print("=" * 80)
    
    # Находим все CSV файлы
    csv_files = sorted(input_dir.glob("merged_ligands_docking_*.csv"))
    
    if not csv_files:
        print(f"❌ Не найдено CSV файлов в {input_dir}")
        return
    
    print(f"\nНайдено {len(csv_files)} файлов для обработки")
    
    for csv_file in csv_files:
        protein = csv_file.stem.replace("merged_ligands_docking_", "").lower()
        print(f"\nОбработка {protein.upper()}...")
        
        try:
            # Читаем CSV
            df = pd.read_csv(csv_file)
            print(f"  Загружено {len(df)} записей")
            
            # Проверяем наличие standard_value
            if "standard_value" not in df.columns:
                print(f"  ⚠️  Отсутствует столбец standard_value, пропускаем")
                continue

            # Фильтруем валидные значения
            valid_mask = df["standard_value"].notna() & (df["standard_value"] > 0)

            # Вычисляем pValue
            # pValue = -log10(Ki в M), где Ki в M = standard_value * 1e-9
            # pValue = -log10(standard_value * 1e-9) = -log10(standard_value) + 9
            # Или: pValue = 9 - log10(standard_value)
            df["pValue"] = np.nan
            df.loc[valid_mask, "pValue"] = 9 - np.log10(df.loc[valid_mask, "standard_value"])
            df["pKi"] = df["pValue"]
            df["is_active"] = valid_mask & (df["standard_value"] < 1000.0)

            # Добавляем класс активности по standard_value (Ki в нМ)
            df["activity_class"] = df["standard_value"].apply(assign_activity_class)

            # Проверяем результат
            n_valid = df["pValue"].notna().sum()
            print(f"  Вычислено pValue для {n_valid} записей")

            if n_valid > 0:
                pvalue_stats = df["pValue"].describe()
                print(
                    f"  pValue: min={pvalue_stats['min']:.2f}, "
                    f"max={pvalue_stats['max']:.2f}, "
                    f"mean={pvalue_stats['mean']:.2f}"
                )

            # Немного статистики по классам активности
            if "activity_class" in df.columns:
                class_counts = df["activity_class"].value_counts(dropna=True).to_dict()
                print(f"  Классы активности (activity_class): {class_counts}")
            
            # Сохраняем в новую директорию
            output_file = output_dir / csv_file.name
            df.to_csv(output_file, index=False)
            print(f"  💾 Сохранен: {output_file.name}")
            
        except Exception as e:
            print(f"  ❌ Ошибка при обработке {csv_file}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'=' * 80}")
    print("ОБРАБОТКА ЗАВЕРШЕНА")
    print(f"{'=' * 80}")
    print(f"📁 Результаты сохранены в: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Добавляет столбец pValue к merged_ligands_docking_*.csv файлам"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help="Директория с исходными CSV файлами (по умолчанию: analysis/tables)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Директория для сохранения файлов с pValue (по умолчанию: analysis/tables)",
    )
    
    args = parser.parse_args()
    
    # Определяем базовую директорию
    base_dir = Path(__file__).resolve().parent.parent.parent
    
    if args.input_dir:
        input_dir = Path(args.input_dir)
    else:
        input_dir = base_dir / "analysis" / "excluding_2z5x_3mjg" / "tables"
    
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = base_dir / "analysis" / "excluding_2z5x_3mjg" / "tables"
    
    # Создаем выходную директорию
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Добавляем столбец pValue
    add_pvalue_column(input_dir, output_dir)


if __name__ == "__main__":
    main()


